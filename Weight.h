#ifndef WEIGHT_H
#define WEIGHT_H

class Weight
{
        public:
                Weight(int index1, int index2, double value) { val = value; i = index1; j = index2; }
                Weight() { Weight(-1, -1, 1.0); }

                void setIndices(int index1, int index2) { i = index1, j = index2;  }
                void setValue(double value) { val = value; }

                int index1() { return i; }
                int index2() { return j; }

                double value() { return val; }

                int operator==(const Weight& aWeight);

        private:
                double val;
                int i;
                int j;

};

int Weight::operator==(const Weight& aWeight)
{
        if ( (i==aWeight.i && j==aWeight.j) || (i==aWeight.j && j==aWeight.i) )
                return 1;
        else
                return 0;
}

#endif
