#ifndef NEIGHBOR_H
#define NEIGHBOR_H

class Neighbor
{
        public:
                Neighbor(int neighborIndex, int weightIndex) { nIndex = neighborIndex; wIndex = weightIndex; }
                Neighbor() { Neighbor(-1, -1); }

                void setNeighborIndex(int neighborIndex) { nIndex = neighborIndex; }
                void setWeightIndex(int weightIndex) { wIndex = weightIndex; }

                int neighborIndex() { return nIndex; }
                int weightIndex() { return wIndex; }

        private:
                int nIndex;
                int wIndex;

};

#endif
