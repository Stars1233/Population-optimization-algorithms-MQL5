//+——————————————————————————————————————————————————————————————————+
//|                                               C_AO_DOAm_dynastic |
//|                                  Copyright 2007-2026, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

#include "#C_AO.mqh"

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
struct S_DOA_Ind
  {
   double c          [];   // координаты (валидные, прогнаны через SeInDiSp)
   double            f;      // фитнес

   void              Init(int coords)
     {
      ArrayResize(c, coords);
      f = -DBL_MAX;
     }
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
class C_AO_DOAm_dynastic : public C_AO
  {
public:
                    ~C_AO_DOAm_dynastic() {}
                     C_AO_DOAm_dynastic()
     {
      ao_name = "DOAm(Dynastic)";
      ao_desc = "Dynastic Optimization Algorithm";
      ao_link = "https://www.mql5.com/ru/articles/23504";

      popSize     = 100;    // размер популяции
      rulerRatio  = 0.05;   // доля правителей rr
      workerRatio = 0.55;   // доля работников rw (исследователи re = 1-rr-rw)
      radW        = 0.4;    // стартовый радиус работников, доля диапазона
      radWmin     = 0.4;    // конечный радиус (== radW -> постоянный, каноника)
      distrPower  = 0.0;    // степень концентрации у правителя (0 -> равномерно, каноника)

      ArrayResize(params, 6);
      params [0].name = "popSize";
      params [0].val = popSize;
      params [1].name = "rulerRatio";
      params [1].val = rulerRatio;
      params [2].name = "workerRatio";
      params [2].val = workerRatio;
      params [3].name = "radW";
      params [3].val = radW;
      params [4].name = "radWmin";
      params [4].val = radWmin;
      params [5].name = "distrPower";
      params [5].val = distrPower;
     }

   void               SetParams()
     {
      popSize     = (int)params [0].val;
      rulerRatio  =      params [1].val;
      workerRatio =      params [2].val;
      radW        =      params [3].val;
      radWmin     =      params [4].val;
      distrPower  =      params [5].val;

      //--- предохранители
      if(popSize < 3)
         popSize = 3;                       // минимум по одному на касту

      if(rulerRatio <= 0.0)
         rulerRatio = 0.01;
      if(rulerRatio > 0.5)
         rulerRatio = 0.5;

      if(workerRatio < 0.0)
         workerRatio = 0.0;
      if(rulerRatio + workerRatio >= 1.0)
         workerRatio = 1.0 - rulerRatio - 0.01;   // исследователям должно что-то остаться
      if(workerRatio < 0.0)
         workerRatio = 0.0;

      if(radW <= 0.0)
         radW = 0.01;
      if(radW > 1.0)
         radW = 1.0;

      if(radWmin <= 0.0)
         radWmin = 0.0001;
      if(radWmin > radW)
         radWmin = radW;

      if(distrPower < 0.0)
         distrPower = 0.0;
     }

   bool               Init(const double &rangeMinP  [],
                           const double &rangeMaxP  [],
                           const double &rangeStepP [],
                           const int     epochsP = 0);

   void               Moving();
   void               Revision();

   //--- видимые параметры
   double             rulerRatio;    // rr
   double             workerRatio;   // rw
   double             radW;          // стартовый радиус работников (доля диапазона)
   double             radWmin;       // конечный радиус (линейное затухание по эпохам)
   double             distrPower;    // степень PowerDistribution (0 -> равномерный куб)

private:
   int                nRulers;       // число правителей
   int                nWorkers;      // число «рядовых» работников (без наследников)
   int                nExplorers;    // число исследователей

   int                epochs;        // всего эпох (от стенда)
   int                epochNow;      // текущая эпоха
   double             radNow;        // актуальный радиус этой эпохи

   S_DOA_Ind          ruler [];      // [nRulers]            — двор: позиции и фитнес
   S_DOA_Ind          pool  [];      // [popSize + nRulers]   — пул для селекции трона
   int                ord   [];      // [popSize + nRulers]   — индексы для сортировки

   //--- вспомогательные
   void               MakeWorker(int slot, int rulerIdx);
   void               MakeExplorer(int slot);
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                              Init                                |
//+------------------------------------------------------------------+
bool C_AO_DOAm_dynastic::Init(const double &rangeMinP  [],
                              const double &rangeMaxP  [],
                              const double &rangeStepP [],
                              const int     epochsP = 0)
  {
   if(!StandardInit(rangeMinP, rangeMaxP, rangeStepP))
      return false;

//--- расчёт каст: rr и rw от popSize, остаток — исследователи
   nRulers = (int)MathRound(rulerRatio * popSize);
   if(nRulers < 1)
      nRulers = 1;
   if(nRulers > popSize - 2)
      nRulers = popSize - 2;

   nWorkers = (int)MathRound(workerRatio * popSize);
   if(nWorkers < 0)
      nWorkers = 0;
   if(nRulers + nWorkers > popSize - 1)
      nWorkers = popSize - 1 - nRulers;    // хотя бы один исследователь

   nExplorers = popSize - nRulers - nWorkers;

//--- расписание радиуса
   epochs   = epochsP;
   if(epochs < 1)
      epochs = 1;
   epochNow = 0;
   radNow   = radW;

//--- буферы (массивы структур)
   ArrayResize(ruler, nRulers);
   ArrayResize(pool,  popSize + nRulers);
   ArrayResize(ord,   popSize + nRulers);

   for(int i = 0; i < nRulers;           i++)
      ruler [i].Init(coords);
   for(int i = 0; i < popSize + nRulers; i++)
      pool  [i].Init(coords);

   return true;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   MakeWorker — кандидат в окрестности radNow вокруг правителя.   |
//|   distrPower = 0 -> равномерный куб (каноника);                  |
//|   distrPower > 0 -> степенная концентрация у правителя (DOAm)    |
//+------------------------------------------------------------------+
void C_AO_DOAm_dynastic::MakeWorker(int slot, int rulerIdx)
  {
   double x, r;

   for(int c = 0; c < coords; c++)
     {
      r = radNow * (rangeMax [c] - rangeMin [c]);

      if(distrPower <= 0.0)
         x = ruler [rulerIdx].c [c] + u.RNDfromCI(-1.0, 1.0) * r;
      else
         x = u.PowerDistribution(ruler [rulerIdx].c [c],
                                 ruler [rulerIdx].c [c] - r,
                                 ruler [rulerIdx].c [c] + r,
                                 distrPower);

      a [slot].c [c] = u.SeInDiSp(x, rangeMin [c], rangeMax [c], rangeStep [c]);
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   MakeExplorer — кандидат равномерно по всему пространству       |
//+------------------------------------------------------------------+
void C_AO_DOAm_dynastic::MakeExplorer(int slot)
  {
   for(int c = 0; c < coords; c++)
      a [slot].c [c] = u.SeInDiSp(u.RNDfromCI(rangeMin [c], rangeMax [c]),
                                  rangeMin [c], rangeMax [c], rangeStep [c]);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                            Moving                                |
//|   Слоты a[i].c — только свежие кандидаты (работники и исследо-   |
//|   ватели); правители живут в ruler[] и стендом не переоцениваются|
//+------------------------------------------------------------------+
void C_AO_DOAm_dynastic::Moving()
  {
//--- первый прогон: вся популяция случайна
   if(!revision)
     {
      for(int i = 0; i < popSize; i++)
         MakeExplorer(i);
      return;
     }

//--- линейное затухание радиуса от radW к radWmin по эпохам
//    (radWmin == radW -> радиус постоянный, каноника)
   epochNow++;
   if(epochNow > epochs)
      epochNow = epochs;

   radNow = radW + (radWmin - radW) * ((double)epochNow / (double)epochs);

   int slot = 0;

//--- 1) наследники: по одному работнику вокруг «своего» правителя
//       (заполняют слоты, освобождённые неоцениваемыми правителями)
   for(int r = 0; r < nRulers; r++)
     {
      MakeWorker(slot, r);
      slot++;
     }

//--- 2) работники: равномерная привязка к правителям по кругу
   for(int j = 0; j < nWorkers; j++)
     {
      MakeWorker(slot, j % nRulers);
      slot++;
     }

//--- 3) исследователи: случайный рестарт в незанятом пространстве
   for(; slot < popSize; slot++)
      MakeExplorer(slot);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                           Revision                               |
//|   Селекция трона: объединённый пул «старые правители + свежие    |
//|   кандидаты», верхние nRulers занимают двор. Правитель уступает  |
//|   место только тому, кто его превзошёл (элитизм по построению).  |
//+------------------------------------------------------------------+
void C_AO_DOAm_dynastic::Revision()
  {
//--- глобальный лучший («император»)
   for(int i = 0; i < popSize; i++)
     {
      if(a [i].f > fB)
        {
         fB = a [i].f;
         ArrayCopy(cB, a [i].c, 0, 0, coords);
        }
     }

//--- пул селекции: свежие кандидаты...
   for(int i = 0; i < popSize; i++)
     {
      ArrayCopy(pool [i].c, a [i].c, 0, 0, coords);
      pool [i].f = a [i].f;
     }

//--- ...плюс действующие правители (на первом проходе их f=-DBL_MAX,
//    и они честно проигрывают любому оценённому кандидату)
   for(int r = 0; r < nRulers; r++)
     {
      ArrayCopy(pool [popSize + r].c, ruler [r].c, 0, 0, coords);
      pool [popSize + r].f = ruler [r].f;
     }

//--- частичная сортировка выбором: нужны только верхние nRulers
   int poolSize = popSize + nRulers;

   for(int i = 0; i < poolSize; i++)
      ord [i] = i;

   for(int i = 0; i < nRulers; i++)
     {
      int bi = i;
      for(int j = i + 1; j < poolSize; j++)
         if(pool [ord [j]].f > pool [ord [bi]].f)
            bi = j;
      int t = ord [i];
      ord [i] = ord [bi];
      ord [bi] = t;
     }

//--- новый двор
   for(int r = 0; r < nRulers; r++)
     {
      ArrayCopy(ruler [r].c, pool [ord [r]].c, 0, 0, coords);
      ruler [r].f = pool [ord [r]].f;
     }

   revision = true;
  }
//+------------------------------------------------------------------+