//+——————————————————————————————————————————————————————————————————+
//|                                                        C_AO_CCEm |
//|                                  Copyright 2007-2026, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

#include "#C_AO.mqh"


//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
struct S_CCEm_Member
  {
   double            c [];
   double            f;

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
class C_AO_CCEm : public C_AO
  {
public:
                    ~C_AO_CCEm() {}
                     C_AO_CCEm()
     {
      ao_name = "CCEm";
      ao_desc = "City Councils Evolution M";
      ao_link = "https://www.mql5.com/ru/articles/24072";

      popSize     = 4;    // всегда 4: кандидаты Y1..Y4 одного члена совета
      councilSize = 20;   // число членов совета (N в статье)
      treeDegree  = 3;    // d - степень дерева советов
      moveProb    = 0.4;  // вероятность участия координаты в ходе

      ArrayResize(params, 4);
      params [0].name = "popSize";
      params [0].val  = popSize;
      params [1].name = "councilSize";
      params [1].val  = councilSize;
      params [2].name = "treeDegree";
      params [2].val  = treeDegree;
      params [3].name = "moveProb";
      params [3].val  = moveProb;
     }

   void              SetParams()
     {
      popSize     = (int)params [0].val;
      councilSize = (int)params [1].val;
      treeDegree  = (int)params [2].val;
      moveProb    =      params [3].val;

      if(councilSize < 2)
         councilSize = 2;
      if(treeDegree < 2)
         treeDegree = 2;

      if(moveProb < 0.0)
         moveProb = 0.0;
      if(moveProb > 1.0)
         moveProb = 1.0;

      params [0].val = popSize;
      params [1].val = councilSize;
      params [2].val = treeDegree;
      params [3].val = moveProb;
     }

   bool              Init(const double &rangeMinP  [],
                          const double &rangeMaxP  [],
                          const double &rangeStepP [],
                          const int     epochsP = 0);

   void              Moving();
   void              Revision();

   //--- видимые параметры
   int               councilSize;
   int               treeDegree;
   double            moveProb;      // вероятность участия координаты в ходе

   //--- диагностика: доля координат кандидатов, упёршихся в границу
   double            ClampRate() { return (probesTotal > 0 ? (double)probesClamped / probesTotal : 0.0); }

private:
   S_CCEm_Member     M  [];
   S_CCEm_Member     Mt [];
   double            partner [];
   bool              gate    [];   // маска двигаемых координат
   int               curMember;
   int               filled;

   long              probesClamped;
   long              probesTotal;

   double            Alpha();
   bool              SameVector(const double &v1 [], const double &v2 []);
   void              MakeGate();
   void              MakeRandomVector(double &v []);
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                              Init                                |
//+------------------------------------------------------------------+
bool C_AO_CCEm::Init(const double &rangeMinP  [],
                     const double &rangeMaxP  [],
                     const double &rangeStepP [],
                     const int     epochsP = 0)
  {
   if(!StandardInit(rangeMinP, rangeMaxP, rangeStepP))
      return false;

   ArrayResize(M,  councilSize);
   ArrayResize(Mt, councilSize);
   ArrayResize(partner, coords);
   ArrayResize(gate,    coords);

   for(int i = 0; i < councilSize; i++)
     {
      M  [i].Init(coords);
      Mt [i].Init(coords);
     }

   curMember     = 0;
   filled        = 0;
   probesClamped = 0;
   probesTotal   = 0;

   return true;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
double C_AO_CCEm::Alpha()
  {
   double rnd = u.RNDprobab();
   return (rnd < 0.5 ? rnd : -rnd);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
bool C_AO_CCEm::SameVector(const double &v1 [], const double &v2 [])
  {
   for(int c = 0; c < coords; c++)
      if(v1 [c] != v2 [c])
         return false;
   return true;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
void C_AO_CCEm::MakeRandomVector(double &v [])
  {
   for(int c = 0; c < coords; c++)
      v [c] = u.SeInDiSp(u.RNDfromCI(rangeMin [c], rangeMax [c]),
                         rangeMin [c], rangeMax [c], rangeStep [c]);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   MakeGate — маска двигаемых координат, общая для Y1..Y4         |
//+------------------------------------------------------------------+
void C_AO_CCEm::MakeGate()
  {
   int cnt = 0;

   for(int c = 0; c < coords; c++)
     {
      gate [c] = (u.RNDprobab() < moveProb);
      if(gate [c])
         cnt++;
     }

//--- ни одной координаты не выбрано - берём одну случайную,
//--- иначе четвёрка выродится в копию текущего решения
   if(cnt == 0)
      gate [u.RNDminusOne(coords)] = true;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                            Moving                                |
//|   Одна эпоха = один член совета.                                 |
//+------------------------------------------------------------------+
void C_AO_CCEm::Moving()
  {
//--- фаза инициализации: заполнение совета случайными решениями
   if(filled < councilSize)
     {
      for(int i = 0; i < popSize; i++)
         MakeRandomVector(a [i].c);
      return;
     }

   int m = curMember;
   int p;

//--- выбор партнёра X1
   if(m == 0)
     {
      //--- корень: партнёр — свежий случайный вектор
      MakeRandomVector(partner);
     }
   else
     {
      //--- родитель узла m в d-арном дереве. Гаусс-Зейдель: p < m всегда.
      p = (m - 1) / treeDegree;

      if(SameVector(M [p].c, M [m].c))
         MakeRandomVector(partner);          // нулевая разность недопустима
      else
         ArrayCopy(partner, M [p].c, 0, 0, coords);
     }

   MakeGate();

//--- четыре кандидата, alpha своя на каждую двигаемую координату
   double x1, x2, alpha, dif, raw;

   for(int c = 0; c < coords; c++)
     {
      if(!gate [c])
        {
         //--- координата не участвует в ходе
         for(int k = 0; k < 4; k++)
            a [k].c [c] = M [m].c [c];
         continue;
        }

      x1    = partner [c];      // X1 - родитель (или случайный партнёр)
      x2    = M [m].c [c];      // X2 - сам член совета
      alpha = Alpha();
      dif   = x1 - x2;

      for(int k = 0; k < 4; k++)
        {
         switch(k)
           {
            case 0:
               raw = x2 + alpha         * dif;
               break;  // к x1
            case 1:
               raw = x2 + (1.0 + alpha) * dif;
               break;  // за x1
            case 2:
               raw = x1 - alpha         * dif;
               break;  // к x2
            default:
               raw = x1 - (1.0 + alpha) * dif;
               break;  // за x2
           }

         probesTotal++;
         if(raw < rangeMin [c] || raw > rangeMax [c])
            probesClamped++;

         a [k].c [c] = u.SeInDiSp(raw, rangeMin [c], rangeMax [c], rangeStep [c]);
        }
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                           Revision                               |
//+------------------------------------------------------------------+
void C_AO_CCEm::Revision()
  {
//--- фаза инициализации
   if(filled < councilSize)
     {
      int order [];
      ArrayResize(order, popSize);
      for(int i = 0; i < popSize; i++)
         order [i] = i;

      for(int i = 0; i < popSize - 1; i++)
         for(int j = i + 1; j < popSize; j++)
            if(a [order [j]].f > a [order [i]].f)
              {
               int t = order [i];
               order [i] = order [j];
               order [j] = t;
              }

      for(int i = 0; i < popSize && filled < councilSize; i++, filled++)
        {
         M [filled].f = a [order [i]].f;
         ArrayCopy(M [filled].c, a [order [i]].c, 0, 0, coords);
        }

      if(filled >= councilSize)
        {
         u.Sorting(M, Mt, councilSize);
         curMember = 0;
        }

      if(M [0].f > fB)
        {
         fB = M [0].f;
         ArrayCopy(cB, M [0].c, 0, 0, coords);
        }

      revision = true;
      return;
     }

//--- лучший из четвёрки, жадный приём
   int    bestIdx = 0;
   double bestF   = a [0].f;

   for(int k = 1; k < 4; k++)
     {
      if(a [k].f > bestF)
        {
         bestF   = a [k].f;
         bestIdx = k;
        }
     }

   int m = curMember;

   if(bestF > M [m].f)
     {
      M [m].f = bestF;
      ArrayCopy(M [m].c, a [bestIdx].c, 0, 0, coords);
     }

   if(M [m].f > fB)
     {
      fB = M [m].f;
      ArrayCopy(cB, M [m].c, 0, 0, coords);
     }

//--- следующий член совета; после полного обхода — перестройка дерева
   curMember++;
   if(curMember >= councilSize)
     {
      curMember = 0;
      u.Sorting(M, Mt, councilSize);
     }
  }
//+------------------------------------------------------------------+
