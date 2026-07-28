//+——————————————————————————————————————————————————————————————————+
//|                                                       C_AO_COBRA |
//|                                  Copyright 2007-2026, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

#include "#C_AO.mqh"

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
struct S_COBRA_Ext                 // дополнительное состояние агента
  {
   double v          [];           // скорость (PSO, BA)
   double            w;            // вес рыбы (FSS)
   int               comp;         // принадлежность компоненту

   void              Init(int coords)
     {
      ArrayResize(v, coords);
      ArrayInitialize(v, 0.0);
      w    = 0.0;
      comp = 0;
     }
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
struct S_COBRA_Comp                // состояние компонента-алгоритма
  {
   double cBest      [];           // лучшие найденные компонентом координаты
   double            fBest;        // лучший фитнес компонента
   double            fAvg;         // средний текущий фитнес субпопуляции
   int               n;            // размер субпопуляции

   void              Init(int coords)
     {
      ArrayResize(cBest, coords);
      ArrayInitialize(cBest, 0.0);
      fBest = -DBL_MAX;
      fAvg  = -DBL_MAX;
      n     = 0;
     }
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
class C_AO_COBRA : public C_AO
  {
public:
                    ~C_AO_COBRA() {}
                     C_AO_COBRA()
     {
      ao_name = "COBRA";
      ao_desc = "Co-Operation of Biology Related Algorithms";
      ao_link = "https://www.mql5.com/ru/articles/23630";

      popSize     = 60;     // суммарный размер популяции (фиксирован)
      minPop      = 2;      // минимальный размер субпопуляции компонента
      commPeriod  = 10;     // период коммуникации, эпох
      commPart    = 0.2;    // доля худших особей, замещаемых при коммуникации
      stagnPeriod = 15;     // эпох застоя до шага ребалансировки

      ArrayResize(params, 5);
      params [0].name = "popSize";
      params [0].val = popSize;
      params [1].name = "minPop";
      params [1].val = minPop;
      params [2].name = "commPeriod";
      params [2].val = commPeriod;
      params [3].name = "commPart";
      params [3].val = commPart;
      params [4].name = "stagnPeriod";
      params [4].val = stagnPeriod;

      //--- константы компонентов (зашиты, см. шапку)
      psoW       = 0.729;   // PSO: инерция (Clerc)
      psoC1      = 1.49445; // PSO: когнитивный коэффициент
      psoC2      = 1.49445; // PSO: социальный коэффициент

      wpsStep    = 1.0;     // WPS: шаг притяжения к вожаку
      wpsBeta    = 0.1;     // WPS: радиус локального поиска вожака (доля диапазона)

      ffaBeta0   = 1.0;     // FFA: базовая привлекательность
      ffaGamma   = 4.0;     // FFA: поглощение (на нормированной дистанции)
      ffaAlpha   = 0.05;    // FFA: случайная составляющая (доля диапазона)

      csaAlpha   = 0.2;     // CSA: масштаб шага Леви
      csaPa      = 0.25;    // CSA: вероятность abandonment (покоординатно)

      baFmin     = 0.0;     // BA: нижняя частота
      baFmax     = 1.0;     // BA: верхняя частота
      baSigma    = 0.01;    // BA: радиус локального поиска (доля диапазона)

      fssStepInd = 0.05;    // FSS: индивидуальный шаг (доля диапазона)
      fssStepVol = 0.1;     // FSS: волевой шаг
      fssWmax    = 10.0;    // FSS: максимальный вес рыбы
     }

   void               SetParams()
     {
      popSize     = (int)params [0].val;
      minPop      = (int)params [1].val;
      commPeriod  = (int)params [2].val;
      commPart    =      params [3].val;
      stagnPeriod = (int)params [4].val;

      //--- предохранители
      if(popSize < compNum)
         popSize = compNum;                 // минимум по одному агенту на компонент

      if(minPop < 1)
         minPop = 1;
      if(minPop > popSize / compNum)
         minPop = popSize / compNum;

      if(commPeriod < 1)
         commPeriod = 1;

      if(commPart < 0.0)
         commPart = 0.0;
      if(commPart > 0.9)
         commPart = 0.9;

      if(stagnPeriod < 1)
         stagnPeriod = 1;
     }

   bool               Init(const double &rangeMinP  [],
                           const double &rangeMaxP  [],
                           const double &rangeStepP [],
                           const int     epochsP = 0);

   void               Moving();
   void               Revision();

   //--- видимые параметры
   int                minPop;        // минимальный размер субпопуляции
   int                commPeriod;    // период коммуникации (эпох)
   double             commPart;      // доля замещаемых худших особей
   int                stagnPeriod;   // эпох застоя до ребалансировки

private:
   static const int   compNum;       // число компонентов (6)

   int                epochs;        // всего эпох (от стенда)
   int                epochNow;      // текущая эпоха

   S_COBRA_Ext        ext  [];       // [popSize] — доп. состояние агентов
   S_COBRA_Comp       comp [];       // [compNum] — состояние компонентов

   //--- агрегаты FSS (по данным предыдущей эпохи)
   double             fssInst [];    // [coords] — инстинктивное смещение косяка
   double             fssBary [];    // [coords] — барицентр косяка
   double             fssWsumPr;     // суммарный вес косяка на предыдущей эпохе
   double             fssVolDir;     // направление волевого оператора: +1 сжатие, -1 расширение

   double             fBprev;        // fB предыдущей эпохи (трекинг застоя)
   int                stagnCnt;      // эпох без улучшения fB

   //--- константы компонентов
   double             psoW, psoC1, psoC2;
   double             wpsStep, wpsBeta;
   double             ffaBeta0, ffaGamma, ffaAlpha;
   double             csaAlpha, csaPa;
   double             baFmin, baFmax, baSigma;
   double             fssStepInd, fssStepVol, fssWmax;

   //--- вспомогательные
   void               MoveAgent(int i, double scale);
   void               Communication();
   void               MigrateAgent(int agentInd, int compTo);
   int                WorstInComp(int k);
   int                RandInComp(int k);
   double             GaussBM();
   double             LevyStep();
  };

const int C_AO_COBRA::compNum = 6;
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                              Init                                |
//+------------------------------------------------------------------+
bool C_AO_COBRA::Init(const double &rangeMinP  [],
                      const double &rangeMaxP  [],
                      const double &rangeStepP [],
                      const int     epochsP = 0)
  {
   if(!StandardInit(rangeMinP, rangeMaxP, rangeStepP))
      return false;

   epochs   = epochsP;
   if(epochs < 1)
      epochs = 1;
   epochNow = 0;

   fBprev   = -DBL_MAX;
   stagnCnt = 0;

//--- буферы (массивы структур)
   ArrayResize(ext,  popSize);
   ArrayResize(comp, compNum);

   for(int i = 0; i < popSize; i++)
      ext [i].Init(coords);
   for(int k = 0; k < compNum; k++)
      comp [k].Init(coords);

   ArrayResize(fssInst, coords);
   ArrayResize(fssBary, coords);
   ArrayInitialize(fssInst, 0.0);
   ArrayInitialize(fssBary, 0.0);

   fssWsumPr = 0.0;
   fssVolDir = 1.0;

//--- равномерное распределение агентов по компонентам
   for(int i = 0; i < popSize; i++)
     {
      ext [i].comp = i % compNum;
      ext [i].w    = fssWmax * 0.5;
      comp [ext [i].comp].n++;
     }

   return true;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                            Moving                                |
//+------------------------------------------------------------------+
void C_AO_COBRA::Moving()
  {
//--- первый прогон: вся популяция случайна
   if(!revision)
     {
      for(int i = 0; i < popSize; i++)
        {
         for(int c = 0; c < coords; c++)
            a [i].c [c] = u.SeInDiSp(u.RNDfromCI(rangeMin [c], rangeMax [c]),
                                     rangeMin [c], rangeMax [c], rangeStep [c]);
        }

      revision = true;
      return;
     }

   epochNow++;
   if(epochNow > epochs)
      epochNow = epochs;

//--- отжиг радиусов локальных поисков (WPS, FFA, BA, FSS)
   double scale = 1.0 - 0.9 * (double)epochNow / (double)epochs;

   for(int i = 0; i < popSize; i++)
      MoveAgent(i, scale);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   MoveAgent — один шаг агента по правилам его компонента         |
//+------------------------------------------------------------------+
void C_AO_COBRA::MoveAgent(int i, double scale)
  {
   int    k = ext [i].comp;
   double x, r1, r2;

   switch(k)
     {
      //--- 0: PSO ----------------------------------------------------
      case 0:
        {
         for(int c = 0; c < coords; c++)
           {
            r1 = u.RNDprobab();
            r2 = u.RNDprobab();

            ext [i].v [c] = psoW  * ext [i].v [c] +
                            psoC1 * r1 * (a [i].cB [c]       - a [i].c [c]) +
                            psoC2 * r2 * (comp [k].cBest [c] - a [i].c [c]);

            a [i].c [c] += ext [i].v [c];
           }
         break;
        }

      //--- 1: WPS ----------------------------------------------------
      case 1:
        {
         bool leader = (a [i].fB == comp [k].fBest);

         for(int c = 0; c < coords; c++)
           {
            if(leader)
               a [i].c [c] += u.RNDfromCI(-1.0, 1.0) * wpsBeta * scale * (rangeMax [c] - rangeMin [c]);
            else
               a [i].c [c] += u.RNDprobab() * wpsStep * (comp [k].cBest [c] - a [i].c [c]);
           }
         break;
        }

      //--- 2: FFA ----------------------------------------------------
      case 2:
        {
         for(int j = 0; j < popSize; j++)
           {
            if(ext [j].comp != k || j == i)
               continue;
            if(a [j].fB <= a [i].fB)
               continue;

            //--- нормированная квадратичная дистанция до более яркого
            double d2 = 0.0;
            for(int c = 0; c < coords; c++)
              {
               double dn = (a [j].cB [c] - a [i].c [c]) / (rangeMax [c] - rangeMin [c]);
               d2 += dn * dn;
              }
            d2 /= (double)coords;

            double beta = ffaBeta0 * MathExp(-ffaGamma * d2);

            for(int c = 0; c < coords; c++)
               a [i].c [c] += beta * (a [j].cB [c] - a [i].c [c]) +
                              ffaAlpha * scale * (u.RNDprobab() - 0.5) * (rangeMax [c] - rangeMin [c]);
           }
         break;
        }

      //--- 3: CSA ----------------------------------------------------
      case 3:
        {
         //--- полёт Леви вокруг лучшего решения компонента
         for(int c = 0; c < coords; c++)
            a [i].c [c] += csaAlpha * LevyStep() * (a [i].c [c] - comp [k].cBest [c]);

         //--- abandonment: смещённое случайное блуждание по подмножеству координат
         int p1 = RandInComp(k);
         int p2 = RandInComp(k);

         for(int c = 0; c < coords; c++)
           {
            if(u.RNDprobab() < csaPa)
               a [i].c [c] += u.RNDprobab() * (a [p1].cB [c] - a [p2].cB [c]);
           }
         break;
        }

      //--- 4: BA -----------------------------------------------------
      case 4:
        {
         double freq = baFmin + (baFmax - baFmin) * u.RNDprobab();

         for(int c = 0; c < coords; c++)
           {
            ext [i].v [c] += (a [i].c [c] - comp [k].cBest [c]) * freq;
            a [i].c [c]   += ext [i].v [c];
           }

         //--- локальный поиск вокруг лучшего решения компонента
         if(u.RNDprobab() > 0.5)
           {
            for(int c = 0; c < coords; c++)
               a [i].c [c] = comp [k].cBest [c] + baSigma * scale * GaussBM() * (rangeMax [c] - rangeMin [c]);
           }
         break;
        }

      //--- 5: FSS ----------------------------------------------------
      default:
        {
         for(int c = 0; c < coords; c++)
           {
            //--- индивидуальное движение
            x = a [i].c [c] + u.RNDfromCI(-1.0, 1.0) * fssStepInd * scale * (rangeMax [c] - rangeMin [c]);

            //--- коллективно-инстинктивное движение (данные предыдущей эпохи)
            x += fssInst [c];

            //--- коллективно-волевое движение (данные предыдущей эпохи)
            x -= fssVolDir * u.RNDprobab() * fssStepVol * scale * (x - fssBary [c]);

            a [i].c [c] = x;
           }
         break;
        }
     }

   for(int c = 0; c < coords; c++)
      a [i].c [c] = u.SeInDiSp(a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                           Revision                               |
//+------------------------------------------------------------------+
void C_AO_COBRA::Revision()
  {
//--- личные и глобальный лучшие
   for(int i = 0; i < popSize; i++)
     {
      if(a [i].f > a [i].fB)
        {
         a [i].fB = a [i].f;
         ArrayCopy(a [i].cB, a [i].c, 0, 0, WHOLE_ARRAY);
        }

      if(a [i].f > fB)
        {
         fB = a [i].f;
         ArrayCopy(cB, a [i].c, 0, 0, WHOLE_ARRAY);
        }
     }

//--- статистика компонентов
   for(int k = 0; k < compNum; k++)
     {
      comp [k].fAvg = 0.0;

      for(int i = 0; i < popSize; i++)
        {
         if(ext [i].comp != k)
            continue;

         comp [k].fAvg += a [i].f;

         if(a [i].fB > comp [k].fBest)
           {
            comp [k].fBest = a [i].fB;
            ArrayCopy(comp [k].cBest, a [i].cB, 0, 0, WHOLE_ARRAY);
           }
        }

      if(comp [k].n > 0)
         comp [k].fAvg /= (double)comp [k].n;
      else
         comp [k].fAvg = -DBL_MAX;
     }

//--- трекинг застоя
   if(fB > fBprev)
     {
      fBprev   = fB;
      stagnCnt = 0;
     }
   else
      stagnCnt++;

//--- миграция к победителю: лучший средний фитнес забирает у худшего
   int winner = 0;
   int loser  = 0;

   for(int k = 1; k < compNum; k++)
     {
      if(comp [k].fAvg > comp [winner].fAvg)
         winner = k;
      if(comp [k].fAvg < comp [loser].fAvg)
         loser = k;
     }

   if(winner != loser && comp [loser].n > minPop)
      MigrateAgent(WorstInComp(loser), winner);

//--- застой: шаг ребалансировки к равномерному распределению
   if(stagnCnt >= stagnPeriod)
     {
      int kMax = 0;
      int kMin = 0;

      for(int k = 1; k < compNum; k++)
        {
         if(comp [k].n > comp [kMax].n)
            kMax = k;
         if(comp [k].n < comp [kMin].n)
            kMin = k;
        }

      if(kMax != kMin && comp [kMax].n > comp [kMin].n + 1)
         MigrateAgent(WorstInComp(kMax), kMin);
     }

//--- коммуникация между популяциями
   if(epochNow > 0 && epochNow % commPeriod == 0)
      Communication();

//--- бухгалтерия FSS (использует f против fP предыдущей эпохи)
   double dfMax = 0.0;
   double df;

   for(int i = 0; i < popSize; i++)
     {
      if(ext [i].comp != 5)
         continue;
      if(a [i].fP == -DBL_MAX)
         continue;

      df = a [i].f - a [i].fP;
      if(df > dfMax)
         dfMax = df;
     }

   double sumDf = 0.0;
   double wSum  = 0.0;

   ArrayInitialize(fssInst, 0.0);
   ArrayInitialize(fssBary, 0.0);

   for(int i = 0; i < popSize; i++)
     {
      if(ext [i].comp != 5)
         continue;

      df = (a [i].fP == -DBL_MAX) ? 0.0 : a [i].f - a [i].fP;
      if(df < 0.0)
         df = 0.0;

      //--- кормление
      if(dfMax > 0.0)
        {
         ext [i].w += df / dfMax;
         if(ext [i].w > fssWmax)
            ext [i].w = fssWmax;
         if(ext [i].w < 1.0)
            ext [i].w = 1.0;
        }

      //--- числитель инстинктивного оператора: успешные смещения
      if(df > 0.0)
        {
         for(int c = 0; c < coords; c++)
            fssInst [c] += df * (a [i].c [c] - a [i].cP [c]);
         sumDf += df;
        }

      //--- числитель барицентра
      for(int c = 0; c < coords; c++)
         fssBary [c] += ext [i].w * a [i].c [c];
      wSum += ext [i].w;
     }

   if(sumDf > 0.0)
     {
      for(int c = 0; c < coords; c++)
         fssInst [c] /= sumDf;
     }
   else
      ArrayInitialize(fssInst, 0.0);

   if(wSum > 0.0)
     {
      for(int c = 0; c < coords; c++)
         fssBary [c] /= wSum;
     }

   fssVolDir = (wSum >= fssWsumPr) ? 1.0 : -1.0;   // косяк потяжелел -> сжатие
   fssWsumPr = wSum;

//--- сохранение предыдущего состояния для следующей эпохи
   for(int i = 0; i < popSize; i++)
     {
      a [i].fP = a [i].f;
      ArrayCopy(a [i].cP, a [i].c, 0, 0, WHOLE_ARRAY);
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   Communication — доля худших особей каждой популяции замещается |
//|   лучшими решениями других компонентов (доноры по кругу)         |
//+------------------------------------------------------------------+
void C_AO_COBRA::Communication()
  {
   for(int k = 0; k < compNum; k++)
     {
      int m = (int)MathRound(commPart * comp [k].n);
      if(m < 1)
         continue;
      if(m > comp [k].n - 1)
         m = comp [k].n - 1;               // лучшую особь не трогаем

      int donor = (k + 1) % compNum;

      for(int r = 0; r < m; r++)
        {
         int w = WorstInComp(k);
         if(w < 0)
            break;

         if(donor == k)
            donor = (donor + 1) % compNum; // свой компонент донором не бывает

         ArrayCopy(a [w].c,  comp [donor].cBest, 0, 0, WHOLE_ARRAY);
         ArrayCopy(a [w].cB, comp [donor].cBest, 0, 0, WHOLE_ARRAY);
         a [w].fB = comp [donor].fBest;
         a [w].f  = comp [donor].fBest;

         ArrayInitialize(ext [w].v, 0.0);
         ext [w].w = fssWmax * 0.5;

         donor = (donor + 1) % compNum;
        }
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   MigrateAgent — смена принадлежности агента (слот сохраняется), |
//|   скоростные состояния мигранта сбрасываются                     |
//+------------------------------------------------------------------+
void C_AO_COBRA::MigrateAgent(int agentInd, int compTo)
  {
   if(agentInd < 0)
      return;

   int compFrom = ext [agentInd].comp;
   if(compFrom == compTo)
      return;

   comp [compFrom].n--;
   comp [compTo].n++;
   ext [agentInd].comp = compTo;

   ArrayInitialize(ext [agentInd].v, 0.0);
   ext [agentInd].w = fssWmax * 0.5;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   WorstInComp — индекс худшего (по текущему f) агента компонента |
//+------------------------------------------------------------------+
int C_AO_COBRA::WorstInComp(int k)
  {
   int    ind = -1;
   double f   = DBL_MAX;

   for(int i = 0; i < popSize; i++)
     {
      if(ext [i].comp != k)
         continue;

      if(a [i].f < f)
        {
         f   = a [i].f;
         ind = i;
        }
     }

   return ind;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   RandInComp — случайный агент внутри компонента                 |
//+------------------------------------------------------------------+
int C_AO_COBRA::RandInComp(int k)
  {
   int cnt = 0;

   for(int i = 0; i < popSize; i++)
      if(ext [i].comp == k)
         cnt++;

   if(cnt == 0)
      return 0;

   int target = u.RNDintInRange(0, cnt - 1);

   for(int i = 0; i < popSize; i++)
     {
      if(ext [i].comp != k)
         continue;
      if(target == 0)
         return i;
      target--;
     }

   return 0;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   GaussBM — стандартная нормальная величина, Бокс-Мюллер         |
//+------------------------------------------------------------------+
double C_AO_COBRA::GaussBM()
  {
   double r1 = u.RNDprobab();
   double r2 = u.RNDprobab();

   if(r1 <= 0.0)
      r1 = 1.0e-12;

   return MathSqrt(-2.0 * MathLog(r1)) * MathCos(2.0 * M_PI * r2);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   LevyStep — шаг полёта Леви, алгоритм Мантеньи, beta = 1.5      |
//+------------------------------------------------------------------+
double C_AO_COBRA::LevyStep()
  {
   static const double beta   = 1.5;
   static const double sigmaU = 0.696575;   // предвычислено для beta = 1.5

   double gU = GaussBM() * sigmaU;
   double gV = GaussBM();

   if(MathAbs(gV) < 1.0e-12)
      gV = 1.0e-12;

   double step = gU / MathPow(MathAbs(gV), 1.0 / beta);

//--- защита от взрывных шагов
   if(step >  10.0)
      step =  10.0;
   if(step < -10.0)
      step = -10.0;

   return step;
  }
//+------------------------------------------------------------------+