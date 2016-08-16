#include <cmath>

#define CELL_MAX 1000000
#define MAX_CHILD 2
#define TREE_DEPTH 2
#define MAX_COPY 2
#define MAX_SUBTYPE ((((int)(pow(MAX_CHILD,TREE_DEPTH)))-1)/(MAX_CHILD-1))
#define NONLEAF ((((int)(pow(MAX_CHILD,TREE_DEPTH-1)))-1)/(MAX_CHILD-1))
#define MAX_H ((int)(pow(2,MAX_CHILD)))
