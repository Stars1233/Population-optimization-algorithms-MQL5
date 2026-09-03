//+——————————————————————————————————————————————————————————————————+
//|                                                        C_AO_CSBO |
//|                                  Copyright 2007-2026, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

#include "#C_AO.mqh"


//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
class C_AO_CSBO : public C_AO
  {
public:
                    ~C_AO_CSBO() {}
                     C_AO_CSBO()
     {
      ao_name = "CSBO";
      ao_desc = "Circulatory System Based Optimization";
      ao_link = "https://www.mql5.com/ru/articles/24504";

      popSize     = 10;       // размер популяции
      nrPart      = 0.0;      // доля худших в малом круге (канон: 1/3)
      crossProb   = 0.1;      // вероятность взять координату мутанта в большом круге
      cauchyScale = 0.2;      // масштаб шага Коши в долях диапазона

      ArrayResize(params, 4);
      params [0].name = "popSize";
      params [0].val = popSize;
      params [1].name = "nrPart";
      params [1].val = nrPart;
      params [2].name = "crossProb";
      params [2].val = crossProb;
      params [3].name = "cauchyScale";
      params [3].val = cauchyScale;
     }

   void               SetParams()
     {
      popSize     = (int)params [0].val;
      nrPart      =      params [1].val;
      crossProb   =      params [2].val;
      cauchyScale =      params [3].val;

      //--- предохранители
      if(popSize < 4)            // агент + три различных партнёра
         popSize = 4;

      if(nrPart < 0.0)
         nrPart = 0.0;
      if(nrPart > 1.0)
         nrPart = 1.0;

      if(crossProb < 0.0)
         crossProb = 0.0;
      if(crossProb > 1.0)
         crossProb = 1.0;

      if(cauchyScale < 0.0)
         cauchyScale = 0.0;

      params [0].val = popSize;
      params [1].val = nrPart;
      params [2].val = crossProb;
      params [3].val = cauchyScale;
     }

   bool               Init(const double &rangeMinP  [],
                           const double &rangeMaxP  [],
                           const double &rangeStepP [],
                           const int     epochsP = 0);

   void               Moving();
   void               Revision();

   //--- видимые параметры
   double             nrPart;         // доля популяции в малом круге
   double             crossProb;      // CR большого круга
   double             cauchyScale;    // масштаб шага малого круга

private:
   double             p [];           // [popSize] — скалярный коэффициент p_i (фаза вен)

   int                nr;             // число агентов в малом круге
   int                phase;          // 0 — вены; 1 — большой + малый круг
   int                feCount;        // выполненных вызовов ФФ (аналог `it` канона)

   bool               nCached;        // кэш Бокса-Мюллера
   double             nSpare;

   //--- вспомогательные
   double             RandN();
   double             RandC();
   void               Pick3(int i, int &a1, int &a2, int &a3);
   void               SwapAgents(int i, int j);
   void               SortByFitness();

   void               MakeVein(int i);
   void               MakeSystemic(int i, double fBestPop, double fWorstPop);
   void               MakePulmonary(int i);
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                              Init                                |
//+------------------------------------------------------------------+
bool C_AO_CSBO::Init(const double &rangeMinP  [],
                     const double &rangeMaxP  [],
                     const double &rangeStepP [],
                     const int     epochsP = 0)
  {
   if(!StandardInit(rangeMinP, rangeMaxP, rangeStepP))
      return false;

   nr = (int)MathRound(popSize * nrPart);
   if(nr < 0)
      nr = 0;
   if(nr > popSize - 1)
      nr = popSize - 1;

   phase   = 0;
   feCount = 0;
   nCached = false;
   nSpare  = 0.0;

   ArrayResize(p, popSize);
   ArrayInitialize(p, 0.0);

   return true;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   RandN — N(0,1) по Боксу-Мюллеру, срез на +-4 сигмы, кэш пары   |
//+------------------------------------------------------------------+
double C_AO_CSBO::RandN()
  {
   if(nCached)
     {
      nCached = false;
      return nSpare;
     }

   double u1 = u.RNDprobab();
   double u2 = u.RNDprobab();

   if(u1 < 1.0e-12)
      u1 = 1.0e-12;

   double rad = MathSqrt(-2.0 * MathLog(u1));
   double ang = 2.0 * M_PI * u2;

   double z0 = rad * MathCos(ang);
   double z1 = rad * MathSin(ang);

   if(z1 >  4.0)
      z1 =  4.0;
   if(z1 < -4.0)
      z1 = -4.0;

   nSpare  = z1;
   nCached = true;

   if(z0 >  4.0)
      z0 =  4.0;
   if(z0 < -4.0)
      z0 = -4.0;

   return z0;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   RandC — стандартное Коши: tan(pi*(u-0.5)), хвост не режем      |
//|   (SeInDiSp сам загонит выброс в диапазон)                       |
//+------------------------------------------------------------------+
double C_AO_CSBO::RandC()
  {
   double v = u.RNDprobab();

   if(v < 1.0e-12)
      v = 1.0e-12;
   if(v > 1.0 - 1.0e-12)
      v = 1.0 - 1.0e-12;

   return MathTan(M_PI * (v - 0.5));
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   Pick3 — три попарно различных индекса, все != i (randperm)     |
//+------------------------------------------------------------------+
void C_AO_CSBO::Pick3(int i, int &a1, int &a2, int &a3)
  {
   do
      a1 = u.RNDminusOne(popSize);
   while(a1 == i);

   do
      a2 = u.RNDminusOne(popSize);
   while(a2 == i || a2 == a1);

   do
      a3 = u.RNDminusOne(popSize);
   while(a3 == i || a3 == a1 || a3 == a2);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   SwapAgents — обмен принятых состояний (cB, fB, p)              |
//+------------------------------------------------------------------+
void C_AO_CSBO::SwapAgents(int i, int j)
  {
   double tmpC [];
   ArrayResize(tmpC, coords);

   ArrayCopy(tmpC,      a [i].cB, 0, 0, coords);
   ArrayCopy(a [i].cB, a [j].cB, 0, 0, coords);
   ArrayCopy(a [j].cB, tmpC,      0, 0, coords);

   double tf  = a [i].fB;
   a [i].fB   = a [j].fB;
   a [j].fB   = tf;

   double tp  = p [i];
   p [i]      = p [j];
   p [j]      = tp;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   SortByFitness — по убыванию fB (лучшие первыми), вставками     |
//+------------------------------------------------------------------+
void C_AO_CSBO::SortByFitness()
  {
   for(int i = 1; i < popSize; i++)
     {
      int j = i;
      while(j > 0 && a [j].fB > a [j - 1].fB)
        {
         SwapAgents(j, j - 1);
         j--;
        }
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   MakeVein — движение крови в венах, (8)(9) канона.              |
//|   K2 = +1 если a1 хуже i (уход от худшего), -1 если лучше        |
//|   (движение к лучшему), при равенстве +1 как в коде авторов.     |
//|   K1 — то же для пары (a2, a3). p_i — скаляр на агента.          |
//+------------------------------------------------------------------+
void C_AO_CSBO::MakeVein(int i)
  {
   int a1, a2, a3;
   Pick3(i, a1, a2, a3);

   double K1 = (a [a2].fB < a [a3].fB) ? 1.0 : (a [a2].fB > a [a3].fB) ? -1.0 : 1.0;
   double K2 = (a [a1].fB < a [i].fB)  ? 1.0 : (a [a1].fB > a [i].fB)  ? -1.0 : 1.0;

   double x;

   for(int c = 0; c < coords; c++)
     {
      x = a [i].cB [c]
          + K2 * p [i] * (a [i].cB [c]  - a [a1].cB [c])
          + K1 * p [i] * (a [a3].cB [c] - a [a2].cB [c]);

      a [i].c [c] = u.SeInDiSp(x, rangeMin [c], rangeMax [c], rangeStep [c]);
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   MakeSystemic — большой круг, (12)(13) канона + кроссовер.      |
//|   ratio = 1 для лучшего, 0 для худшего; запоминается в p_i для   |
//|   следующей фазы вен. Мутант x_a1 + ratio*(x_a3 - x_a2), из      |
//|   него берётся координата с вероятностью crossProb.              |
//+------------------------------------------------------------------+
void C_AO_CSBO::MakeSystemic(int i, double fBestPop, double fWorstPop)
  {
   double spread = fBestPop - fWorstPop;
   double ratio  = (spread > 1.0e-12) ? (a [i].fB - fWorstPop) / spread : 1.0;

   p [i] = ratio;

   int a1, a2, a3;
   Pick3(i, a1, a2, a3);

   double x;

   for(int c = 0; c < coords; c++)
     {
      if(u.RNDprobab() < crossProb)
        {
         x = a [a1].cB [c] + ratio * (a [a3].cB [c] - a [a2].cB [c]);
         a [i].c [c] = u.SeInDiSp(x, rangeMin [c], rangeMax [c], rangeStep [c]);
        }
      else
         a [i].c [c] = a [i].cB [c];
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   MakePulmonary — малый круг, (10) канона.                       |
//|   Один скаляр N(0,1)/feCount на агента, вектор Коши по           |
//|   координатам, масштаб — доля ширины диапазона.                  |
//+------------------------------------------------------------------+
void C_AO_CSBO::MakePulmonary(int i)
  {
   double s = cauchyScale * RandN() / (double)MathMax(feCount, 1);
   double rSize, x;

   for(int c = 0; c < coords; c++)
     {
      rSize = rangeMax [c] - rangeMin [c];

      x = a [i].cB [c] + s * RandC() * rSize;

      a [i].c [c] = u.SeInDiSp(x, rangeMin [c], rangeMax [c], rangeStep [c]);
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                            Moving                                |
//+------------------------------------------------------------------+
void C_AO_CSBO::Moving()
  {
//--- первый прогон: стартовая популяция
   if(!revision)
     {
      for(int i = 0; i < popSize; i++)
         for(int c = 0; c < coords; c++)
            a [i].c [c] = u.SeInDiSp(u.RNDfromCI(rangeMin [c], rangeMax [c]),
                                     rangeMin [c], rangeMax [c], rangeStep [c]);
      return;
     }

//--- фаза 0: вены, все агенты
   if(phase == 0)
     {
      for(int i = 0; i < popSize; i++)
         MakeVein(i);
      return;
     }

//--- фаза 1: популяция отсортирована в Revision фазы 0
   double fBestPop  = -DBL_MAX;
   double fWorstPop =  DBL_MAX;

   for(int i = 0; i < popSize; i++)
     {
      if(a [i].fB > fBestPop)
         fBestPop  = a [i].fB;
      if(a [i].fB < fWorstPop)
         fWorstPop = a [i].fB;
     }

   int nl = popSize - nr;

   for(int i = 0; i < nl; i++)
      MakeSystemic(i, fBestPop, fWorstPop);

   for(int i = nl; i < popSize; i++)
      MakePulmonary(i);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                           Revision                               |
//+------------------------------------------------------------------+
void C_AO_CSBO::Revision()
  {
//--- глобальный лучший
   for(int i = 0; i < popSize; i++)
     {
      if(a [i].f > fB)
        {
         fB = a [i].f;
         ArrayCopy(cB, a [i].c, 0, 0, coords);
        }
     }

   feCount += popSize;

//--- первый проход: стартовая популяция становится принятой
   if(!revision)
     {
      for(int i = 0; i < popSize; i++)
        {
         ArrayCopy(a [i].cB, a [i].c, 0, 0, coords);
         a [i].fB = a [i].f;
         p [i]    = u.RNDprobab();
        }

      phase    = 0;
      revision = true;
      return;
     }

//--- жадная приёмка (строго лучше, как в коде авторов)
   for(int i = 0; i < popSize; i++)
     {
      if(a [i].f > a [i].fB)
        {
         a [i].fB = a [i].f;
         ArrayCopy(a [i].cB, a [i].c, 0, 0, coords);
        }
     }

   if(phase == 0)
     {
      //--- после вен: сортировка, дальше большой + малый круг
      SortByFitness();
      phase = 1;
      return;
     }

//--- после малого круга: сброс p_i у худших (11), обратно в вены
   for(int i = popSize - nr; i < popSize; i++)
      p [i] = u.RNDprobab();

   phase = 0;
  }
//+------------------------------------------------------------------+