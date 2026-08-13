//+——————————————————————————————————————————————————————————————————+
//|                                                        C_AO_FDOm |
//|                                  Copyright 2007-2026, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

#include "#C_AO.mqh"

//+------------------------------------------------------------------+
//| Хранит только то, чего нет в S_AO_Agent. Принятая позиция агента |
//| и её фитнес живут в a[i].cB и a[i].fB (при жадном приёме принятая|
//| позиция и есть личный рекорд), состояние автомата — в a[i].cnt.  |
//+------------------------------------------------------------------+ 
struct S_FDOm_Pace
  {
   double            cand [];   // pace текущей основной пробы
   double            prev [];   // сохранённый pace последнего принятого хода
   bool              hasPrev;   // есть ли сохранённая инерция

   void              Init(int coords)
     {
      ArrayResize(cand, coords);
      ArrayResize(prev, coords);
      ArrayInitialize(cand, 0.0);
      ArrayInitialize(prev, 0.0);
      hasPrev = false;
     }
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
class C_AO_FDOm : public C_AO
  {
public:
                    ~C_AO_FDOm() {}
                     C_AO_FDOm()
     {
      ao_name = "FDOm";
      ao_desc = "Fitness Dependent Optimizer M";
      ao_link = "https://www.mql5.com/ru/articles/21639";

      popSize = 50;     // размер популяции (в статье 30 разведчиков)
      wf      = 0.0;    // весовой фактор, каноника: 0 или 1

      ArrayResize(params, 2);
      params [0].name = "popSize";
      params [0].val = popSize;
      params [1].name = "wf";
      params [1].val = wf;
     }

   void               SetParams()
     {
      popSize = (int)params [0].val;
      wf      =      params [1].val;

      if(popSize < 2)
         popSize = 2;

      //--- в каноне wf принимает только 0 или 1
      if(wf != 0.0 && wf != 1.0)
         wf = (wf >= 0.5) ? 1.0 : 0.0;
     }

   bool               Init(const double &rangeMinP  [],
                           const double &rangeMaxP  [],
                           const double &rangeStepP [],
                           const int     epochsP = 0);

   void               Moving();
   void               Revision();

   //--- видимые параметры
   double             wf;            // весовой фактор из ур. (2)/(6)

private:
   S_FDOm_Pace        pace [];       // [popSize] - векторы шага, которых нет в S_AO_Agent

   //--- состояния автомата, хранятся в a[i].cnt
   //    STAGE_MAIN - основная проба по ур. (3)/(4)/(5)
   //    STAGE_PREV - повторная проба сохранённым pace
#define            STAGE_MAIN 0
#define            STAGE_PREV 1
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                              Init                                |
//+------------------------------------------------------------------+
bool C_AO_FDOm::Init(const double &rangeMinP  [],
                     const double &rangeMaxP  [],
                     const double &rangeStepP [],
                     const int     epochsP = 0)
  {
//--- StandardInit создаёт a[], обнуляет a[i].cnt и ставит a[i].fB = -DBL_MAX
   if(!StandardInit(rangeMinP, rangeMaxP, rangeStepP))
      return false;

   ArrayResize(pace, popSize);
   for(int i = 0; i < popSize; i++)
      pace [i].Init(coords);

   return true;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                            Moving                                |
//|   STAGE_MAIN: основная проба, pace по ур. (3)/(4)/(5)            |
//|   STAGE_PREV: повторная проба сохранённым pace (второй шанс)     |
//|                                                                  |
//|   Условная вторая проба оригинала разложена в автомат на два     |
//|   такта: каждый агент отдаёт стенду ровно одного кандидата за    |
//|   эпоху, поэтому бюджет вычислений ЦФ расходуется точно.         |
//+------------------------------------------------------------------+
void C_AO_FDOm::Moving()
  {
//--- первый прогон: вся популяция случайна
   if(!revision)
     {
      for(int i = 0; i < popSize; i++)
         for(int c = 0; c < coords; c++)
            a [i].c [c] = u.SeInDiSp(u.RNDfromCI(rangeMin [c], rangeMax [c]),
                                     rangeMin [c], rangeMax [c], rangeStep [c]);
      return;
     }

   const double eps = 1e-12;
   double fw, r, p, x, xc;

   for(int i = 0; i < popSize; i++)
     {
      if(a [i].cnt == STAGE_MAIN)
        {
         //--- вес приспособленности, ур. (6) — версия для максимизации.
         //    Делитель здесь fB, поэтому защита от деления на ноль
         //    переезжает с фитнеса агента на фитнес лучшего.
         if(MathAbs(fB) < eps)
            fw = 0.0;
         else
            fw = MathAbs(a [i].fB / fB) - wf;

         //--- вырожденность определяется скалярным fw, значит решается
         //    один раз на агента, а не на координату
         bool degenerate = (MathAbs(fw) < eps || MathAbs(fw - 1.0) < eps);

         for(int c = 0; c < coords; c++)
           {
            //--- ключевое отличие: свой случайный скаляр на каждую координату
            r  = u.RNDfromCI(-1.0, 1.0);
            xc = a [i].cB [c];

            if(degenerate)
               p = xc * r;                                                // ур. (3)
            else
              {
               if(r < 0.0)
                  p = -(xc - cB [c]) * fw;                                // ур. (4)
               else
                  p = (xc - cB [c]) * fw;                                 // ур. (5)
              }

            pace [i].cand [c] = p;

            x = xc + p;                                                   // ур. (1)
            a [i].c [c] = u.SeInDiSp(x, rangeMin [c], rangeMax [c], rangeStep [c]);
           }
        }
      else
        {
         //--- «продолжить в прежнем направлении»
         for(int c = 0; c < coords; c++)
           {
            x = a [i].cB [c] + pace [i].prev [c];
            a [i].c [c] = u.SeInDiSp(x, rangeMin [c], rangeMax [c], rangeStep [c]);
           }
        }
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                           Revision                               |
//|   Жадный приём: позиция агента может только улучшаться, поэтому  |
//|   a[i].cB/fB и есть его текущее положение. Персонального лучшего |
//|   в отдельном смысле в FDO нет, память агента — только pace.     |
//+------------------------------------------------------------------+
void C_AO_FDOm::Revision()
  {
//--- глобальный лучший (в статье обновляется асинхронно, здесь пакетно)
   for(int i = 0; i < popSize; i++)
     {
      if(a [i].f > fB)
        {
         fB = a [i].f;
         ArrayCopy(cB, a [i].c, 0, 0, coords);
        }
     }

//--- первый заход: фиксируем стартовые позиции разведчиков
   if(!revision)
     {
      for(int i = 0; i < popSize; i++)
        {
         ArrayCopy(a [i].cB, a [i].c, 0, 0, coords);
         a [i].fB  = a [i].f;
         a [i].cnt = STAGE_MAIN;

         pace [i].hasPrev = false;
        }

      revision = true;
      return;
     }

   for(int i = 0; i < popSize; i++)
     {
      if(a [i].f > a [i].fB)
        {
         //--- ход принят
         ArrayCopy(a [i].cB, a [i].c, 0, 0, coords);
         a [i].fB = a [i].f;

         if(a [i].cnt == STAGE_MAIN)
           {
            //--- новый pace запоминается для потенциального переиспользования
            ArrayCopy(pace [i].prev, pace [i].cand, 0, 0, coords);
            pace [i].hasPrev = true;
           }
         //--- в STAGE_PREV prev не обновляется: инерция переиспользуется как есть

         a [i].cnt = STAGE_MAIN;
        }
      else
        {
         //--- основная проба провалилась -> следующий такт по сохранённой инерции;
         //    если инерции нет или провалилась уже она — агент остаётся на месте
         if(a [i].cnt == STAGE_MAIN && pace [i].hasPrev)
            a [i].cnt = STAGE_PREV;
         else
            a [i].cnt = STAGE_MAIN;
        }
     }
  }
//+------------------------------------------------------------------+