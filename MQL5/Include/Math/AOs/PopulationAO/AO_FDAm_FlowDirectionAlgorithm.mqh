//+——————————————————————————————————————————————————————————————————+
//|                                                        C_AO_FDAm |
//|                                  Copyright 2007-2026, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

#include "#C_AO.mqh"

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
struct S_FDAm_Flow
  {
   double c          [];    // позиция потока (валидная, прогнана через SeInDiSp)
   double bnC        [];    // позиция лучшего соседа текущей итерации
   double            f;     // фитнес потока
   double            bnF;   // фитнес лучшего соседа текущей итерации

   void              Init(int coords)
     {
      ArrayResize(c,   coords);
      ArrayResize(bnC, coords);
      f   = -DBL_MAX;
      bnF = -DBL_MAX;
     }
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
class C_AO_FDAm : public C_AO
  {
public:
                    ~C_AO_FDAm() {}
                     C_AO_FDAm()
     {
      ao_name = "FDAm";
      ao_desc = "Flow Direction Algorithm M";
      ao_link = "https://www.mql5.com/ru/articles/23731";

      popSize    = 25;     // размер популяции (== число потоков)
      nNeighbors = 1;      // соседей на поток
      radMin     = 0.7;    // пол радиуса зонда, доля диапазона
      distrPower = 50.0;   // концентрация зонда у потока (PowerDistribution)
      useSlope   = 0;      // 0 -> ход только веткой (9); 1 -> каноническое ядро (6)(7)(8)
      moveProb   = 0.5;    // вероятность перемещения координаты в ветке (9)

      ArrayResize(params, 6);
      params [0].name = "popSize";
      params [0].val = popSize;
      params [1].name = "nNeighbors";
      params [1].val = nNeighbors;
      params [2].name = "radMin";
      params [2].val = radMin;
      params [3].name = "distrPower";
      params [3].val = distrPower;
      params [4].name = "useSlope";
      params [4].val = useSlope;
      params [5].name = "moveProb";
      params [5].val = moveProb;
     }

   void               SetParams()
     {
      popSize    = (int)params [0].val;
      nNeighbors = (int)params [1].val;
      radMin     =      params [2].val;
      distrPower =      params [3].val;
      useSlope   = (int)params [4].val;
      moveProb   =      params [5].val;

      //--- предохранители
      if(popSize < 2)
         popSize = 2;

      if(nNeighbors < 1)
         nNeighbors = 1;
      if(nNeighbors > 10)
         nNeighbors = 10;

      if(radMin < 0.0001)
         radMin = 0.0001;
      if(radMin > 1.0)
         radMin = 1.0;

      if(distrPower < 1.0)
         distrPower = 1.0;

      if(useSlope != 1)
         useSlope = 0;

      if(moveProb < 0.01)
         moveProb = 0.01;
      if(moveProb > 1.0)
         moveProb = 1.0;
     }

   bool               Init(const double &rangeMinP  [],
                           const double &rangeMaxP  [],
                           const double &rangeStepP [],
                           const int     epochsP = 0);

   void               Moving();
   void               Revision();

   //--- видимые параметры
   int                nNeighbors;     // число соседей на поток
   double             radMin;         // пол радиуса зонда (доля диапазона)
   double             distrPower;     // концентрация PowerDistribution
   int                useSlope;       // 1 -> каноническое ядро хода (для A/B)
   double             moveProb;       // вероятность перемещения координаты в (9)

private:
   S_FDAm_Flow        flow [];        // [popSize] — потоки

   int                phase;          // 0..nNeighbors-1 сосед, nNeighbors — ход
   int                epochs;         // всего эпох (от стенда)
   int                epochNow;       // текущая эпоха
   double             fSpread;        // размах фитнеса стартовой популяции (для useSlope=1)

   bool               nCached;        // кэш второго значения Бокса-Мюллера
   double             nSpare;

   //--- вспомогательные
   double             RandN();
   void               MakeNeighbor(int i);
   void               MakeMove(int i);
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                              Init                                |
//+------------------------------------------------------------------+
bool C_AO_FDAm::Init(const double &rangeMinP  [],
                     const double &rangeMaxP  [],
                     const double &rangeStepP [],
                     const int     epochsP = 0)
  {
   if(!StandardInit(rangeMinP, rangeMaxP, rangeStepP))
      return false;

   epochs = epochsP;
   if(epochs < 1)
      epochs = 1;

   epochNow = 0;
   phase    = 0;
   fSpread  = 1.0;
   nCached  = false;
   nSpare   = 0.0;

   ArrayResize(flow, popSize);
   for(int i = 0; i < popSize; i++)
      flow [i].Init(coords);

   return true;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   RandN — N(0,1) по Боксу-Мюллеру, срез на +-4 сигмы, кэш пары   |
//+------------------------------------------------------------------+
double C_AO_FDAm::RandN()
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
//|   MakeNeighbor — зонд вокруг потока (M2):                        |
//|   радиус = max(distN, radMin) по каждой координате, степенная    |
//|   концентрация у потока. distN — нормированное расстояние до     |
//|   глобального лучшего (родная связь FDA из (3)).                 |
//+------------------------------------------------------------------+
void C_AO_FDAm::MakeNeighbor(int i)
  {
   double rSize, d, x, r;

//--- нормированное расстояние до глобального лучшего
   double distN = 0.0;

   for(int c = 0; c < coords; c++)
     {
      rSize = rangeMax [c] - rangeMin [c];
      if(rSize <= 0.0)
         continue;

      d = (cB [c] - flow [i].c [c]) / rSize;
      distN += d * d;
     }
   distN = MathSqrt(distN);

//--- адаптивный радиус: далёкий поток шарит широко, близкий дожимает;
//    пол radMin не даёт лучшему потоку замереть (у него distN = 0)
   double radN = distN;
   if(radN < radMin)
      radN = radMin;
   if(radN > 1.0)
      radN = 1.0;

   for(int c = 0; c < coords; c++)
     {
      rSize = rangeMax [c] - rangeMin [c];
      r     = radN * rSize;

      x = u.PowerDistribution(flow [i].c [c],
                              flow [i].c [c] - r,
                              flow [i].c [c] + r,
                              distrPower);

      a [i].c [c] = u.SeInDiSp(x, rangeMin [c], rangeMax [c], rangeStep [c]);
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   MakeMove — ход потока.                                         |
//|   useSlope=0 (M1): всегда ветка (9) — сток к случайному более    |
//|   выгодному потоку либо к глобальному лучшему, в разностной      |
//|   форме, покоординатно.                                          |
//|   useSlope=1: каноническое ядро (6)(7)(8) для A/B-замера         |
//|   (нормировки и срезки — как в честной реализации FDA).          |
//+------------------------------------------------------------------+
void C_AO_FDAm::MakeMove(int i)
  {
   double rSize, dXn, x;

//--- каноническое ядро хода — только по запросу
   if(useSlope == 1 && flow [i].bnF > flow [i].f)
     {
      double Ln = 0.0;

      for(int c = 0; c < coords; c++)
        {
         rSize = rangeMax [c] - rangeMin [c];
         if(rSize <= 0.0)
            continue;

         dXn = (flow [i].c [c] - flow [i].bnC [c]) / rSize;
         Ln += dXn * dXn;
        }
      Ln = MathSqrt(Ln);

      if(Ln > 1.0e-10)
        {
         double dFn = (flow [i].bnF - flow [i].f) / fSpread;
         if(dFn > 1.0)
            dFn = 1.0;

         for(int c = 0; c < coords; c++)
           {
            rSize = rangeMax [c] - rangeMin [c];

            dXn = (flow [i].c [c] - flow [i].bnC [c]) / rSize;

            double S0    = dFn / MathMax(MathAbs(dXn), 1.0e-10);
            double V     = RandN() * S0;
            double stepN = V * (dXn / Ln);

            if(stepN >  1.0)
               stepN =  1.0;
            if(stepN < -1.0)
               stepN = -1.0;

            x = flow [i].c [c] + stepN * rSize;

            a [i].c [c] = u.SeInDiSp(x, rangeMin [c], rangeMax [c], rangeStep [c]);
           }
         return;
        }
     }

//--- (9) как единственный ход (M1): сток к случайному более выгодному
//    потоку, иначе — к глобальному лучшему
   int r = u.RNDminusOne(popSize);

//--- Частичное перемещение (M4): каждая координата двигается с
//    вероятностью moveProb, остальные остаются на месте. Полномерный
//    шаг в высокой размерности почти всегда портит больше координат,
//    чем улучшает, и жадная приёмка его отбрасывает — потоки в
//    1000-мерных тестах жили на одних зондах. При moveProb=1.0
//    поведение прежнее бит-в-бит.
   if(flow [r].f > flow [i].f)
     {
      for(int c = 0; c < coords; c++)
        {
         if(moveProb < 1.0 && u.RNDprobab() > moveProb)
           {
            a [i].c [c] = flow [i].c [c];
            continue;
           }

         x = flow [i].c [c] + RandN() * (flow [r].c [c] - flow [i].c [c]);
         a [i].c [c] = u.SeInDiSp(x, rangeMin [c], rangeMax [c], rangeStep [c]);
        }
     }
   else
     {
      for(int c = 0; c < coords; c++)
        {
         if(moveProb < 1.0 && u.RNDprobab() > moveProb)
           {
            a [i].c [c] = flow [i].c [c];
            continue;
           }

         x = flow [i].c [c] + 2.0 * u.RNDprobab() * (cB [c] - flow [i].c [c]);
         a [i].c [c] = u.SeInDiSp(x, rangeMin [c], rangeMax [c], rangeStep [c]);
        }
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                            Moving                                |
//+------------------------------------------------------------------+
void C_AO_FDAm::Moving()
  {
//--- первый прогон: стартовые позиции потоков
   if(!revision)
     {
      for(int i = 0; i < popSize; i++)
         for(int c = 0; c < coords; c++)
            a [i].c [c] = u.SeInDiSp(u.RNDfromCI(rangeMin [c], rangeMax [c]),
                                     rangeMin [c], rangeMax [c], rangeStep [c]);
      return;
     }

   epochNow++;
   if(epochNow > epochs)
      epochNow = epochs;

//--- фаза разведки окрестности
   if(phase < nNeighbors)
     {
      for(int i = 0; i < popSize; i++)
         MakeNeighbor(i);
      return;
     }

//--- фаза перемещения
   for(int i = 0; i < popSize; i++)
      MakeMove(i);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                           Revision                               |
//+------------------------------------------------------------------+
void C_AO_FDAm::Revision()
  {
//--- глобальный лучший — по всем слотам
   for(int i = 0; i < popSize; i++)
     {
      if(a [i].f > fB)
        {
         fB = a [i].f;
         ArrayCopy(cB, a [i].c, 0, 0, coords);
        }
     }

//--- первый проход: стартовая популяция становится потоками
   if(!revision)
     {
      double fMax = -DBL_MAX;
      double fMin =  DBL_MAX;

      for(int i = 0; i < popSize; i++)
        {
         ArrayCopy(flow [i].c, a [i].c, 0, 0, coords);
         flow [i].f   = a [i].f;
         flow [i].bnF = -DBL_MAX;

         if(a [i].f > fMax)
            fMax = a [i].f;
         if(a [i].f < fMin)
            fMin = a [i].f;
        }

      fSpread = fMax - fMin;
      if(fSpread < 1.0e-10)
         fSpread = 1.0e-10;

      phase    = 0;
      revision = true;
      return;
     }

//--- фаза разведки: копим лучшего соседа и СРАЗУ принимаем зонд,
//    если он лучше потока (M3): оценка уже оплачена бюджетом
   if(phase < nNeighbors)
     {
      for(int i = 0; i < popSize; i++)
        {
         if(a [i].f > flow [i].bnF)
           {
            flow [i].bnF = a [i].f;
            ArrayCopy(flow [i].bnC, a [i].c, 0, 0, coords);
           }

         if(a [i].f > flow [i].f)
           {
            flow [i].f = a [i].f;
            ArrayCopy(flow [i].c, a [i].c, 0, 0, coords);
           }
        }

      phase++;
      return;
     }

//--- фаза перемещения: жадная приёмка хода
   for(int i = 0; i < popSize; i++)
     {
      if(a [i].f > flow [i].f)
        {
         flow [i].f = a [i].f;
         ArrayCopy(flow [i].c, a [i].c, 0, 0, coords);
        }

      flow [i].bnF = -DBL_MAX;
     }

   phase = 0;
  }
//+------------------------------------------------------------------+