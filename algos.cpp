#include <random>
#include <vector>
#include <stack>
#include <tuple>

const short dimen=256;

std::vector<std::vector<short>> map;
std::vector<std::vector<short>> marker;

struct Point {
  short x,y;
  Point(short xx, short yy):x(xx),y(yy){}
};

class LSystem{

public:
    enum Direction {L=0,T,R,B};
    struct State{
        double x, y;
        Direction ang;
    };

    std::stack<State> stack;
    Point curr_pos;
    Direction curr_ang;

    State currState(){
        State s;
        s.x = curr_pos.x;
        s.y = curr_pos.y;
        s.ang = curr_ang;
        return s;
    }

    void forward(){
        Direction len = curr_ang;
        if(len == L)
            curr_pos.x--;
        if(len == R)
            curr_pos.x++;
        if(len == B)
            curr_pos.y--;
        if(len == T)
            curr_pos.y++;

    }

    void rotate(Direction ang){
        if(ang == R){
            switch (curr_ang){
            case L:
                curr_ang = T;break;
            case T:
                curr_ang = R;break;
            case R:
                curr_ang = B;break;
            case B:
                curr_ang = L;break;

            }
        }
        else
            switch (curr_ang){
            case L:
                curr_ang = B;break;
            case T:
                curr_ang = L;break;
            case R:
                curr_ang = T;break;
            case B:
                curr_ang = R;break;

            }
    }

    void push(){
        State s;
        s.x = curr_pos.x;
        s.y = curr_pos.y;
        s.ang = curr_ang;
        stack.push(s);
    }

    void pop(){
        State s = stack.top();
        stack.pop();
        curr_pos.x = s.x;
        curr_pos.y = s.y;
        curr_ang = s.ang;
    }

public:
  LSystem():curr_pos(-1,-1){}

    void setCurrState(State s){
        curr_pos.x = s.x;
        curr_pos.y = s.y;
        curr_ang = R;
    }

    std::tuple< std::vector<State>, std::vector<Point> > iterate(int counter){
        std::vector<Point> draw_points;
        std::vector<State> active_points;
        auto func = [&](){
            auto c = counter;
            while(c--){
                forward();
                draw_points.push_back(curr_pos);
            }
        };
        func();
        push();
        rotate(L);
        func();
        active_points.push_back(currState());
        pop();
        push();
        rotate(R);
        func();
        active_points.push_back(currState());
        pop();

        return std::tie(active_points, draw_points);
    }

};


void DoLSystem(){
    std::vector<LSystem::State> P;
    LSystem::State sp;
    sp.x = dimen/2;
    sp.y = dimen/2;
    sp.ang = LSystem::R;
    marker[dimen/2][dimen/2] = 1;
    P.push_back(sp);
    LSystem l;
    int len = 32;
    for(int i =7; i> 1; --i){
        std::vector<LSystem::State> newPs;
        for(auto& n : P){
            l.setCurrState(n);
            std::vector<LSystem::State> newnodes;
            std::vector<Point> drawnodes;
            std::tie( newnodes, drawnodes ) = l.iterate(len);
            for(auto& xx : drawnodes){
                marker[(int)xx.x][(int)xx.y] = 1;
            }
            for(auto& xx : newnodes){
                newPs.push_back(xx);
            }
        }
        P = newPs;
        newPs.clear();
        len /= 2;
    }

}

template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}

class PoissonDiscSampler{
public:
    typedef std::vector<std::vector<int>> IntMap;
    double radius, radius2;
    int maxTrys;
    int dimen;
    IntMap grid;
    std::vector<Point> activeList;
    std::random_device rd;
    std::mt19937 mt;
    std::uniform_real_distribution<double> rand_angle;
    std::uniform_real_distribution<double> rand_dist;

    PoissonDiscSampler(double rad):
        radius(rad), radius2(rad*2),
        maxTrys(30), dimen(256),
        grid(dimen, std::vector<int>(dimen , -1)),
        mt(rd()),rand_angle(0, 3.14*2), rand_dist(rad, rad*2)
    {
      Point p(dimen/2, dimen/2);
      //p.x = dimen/2; p.y = dimen/2;
        activeList.push_back(p);
        grid[p.x][p.y] = 1;
    }

    bool anypoint(const int x, const int y){
        for(int x1 = x - radius*3; x1 < x + radius*3; ++x1){
            for(int y1 = y - radius*3; y1 < y + radius*3; ++y1){
                if(x1 >= 0 && x1 < dimen && y1 >= 0 && y1 < dimen)
                    if(grid[x1][y1] != -1)
                        if((x1 -x)*(x1-x) + (y1-y)*(y1-y) < radius*radius)
                            return true;
            }
        }
        return false;
    }

    void Iterate(){
        while(!activeList.empty()){
            Point p = *select_randomly(activeList.begin(), activeList.end());
            activeList.erase(std::remove_if(activeList.begin(), activeList.end(),
                                            [&](Point t){return t.x == p.x && t.y == p.y;}), activeList.end());
            for(int i =0; i < maxTrys; ++i){
                double a = rand_angle(mt);
                double r = rand_dist(mt);
                int x = p.x + r * std::cos(a);
                int y = p.y + r * std::sin(a);
                Point n(x,y);// n.x = x; n.y = y;

                if(x>=0 && x < dimen && y >= 0 && y < dimen && !anypoint(x,y)){
                    activeList.push_back(n);
                    grid[x][y] = 1;
                    //break;
                }
            }
        }
    }

};

void DoPoission(){
    PoissonDiscSampler rad(8);
    rad.Iterate();


    for(int x = 0; x < dimen; ++x){
        for(int y = 0; y<dimen; ++y){
	  if(rad.grid[x][y] == 1){
                marker[x][y] =1;
            }
        }
    }
}


class SpaceColonization{
public:
    struct Vector{
        double x,y;
        Vector():x(0),y(0){}
    };

    struct Node : Point {
      std::vector<Point> attractors;
      Vector direction;
      Node():Point(-1,-1){}
    };

    struct AttractionPoint : Point {
      int closestNode;
      AttractionPoint():closestNode(-1), Point(-1,-1){}
      AttractionPoint(int x1, int y1):Point(x1,y1), closestNode(-1){}
    };

    void AddAttractionPoint(int x, int y){
        attractionPoints.push_back( AttractionPoint(x,y));
    }

    void AddInitalNode(int x, int y){
        Node n;
        n.x = x; n.y = y;
        nodes.push_back(n);
    }

    std::vector<Node>&  GetNodes(){
        return nodes;
    }

    double distanceSquare(Point p1, Point p2){
        return (p1.x - p2.x)*(p1.x - p2.x) + (p1.y-p2.y)*(p1.y-p2.y);
    }

    Vector normalizeVec(Vector v){
        double n = std::sqrt(v.x * v.x + v.y * v.y);
        v.x = v.x/n;
        v.y = v.y/n;
        return v;
    }

    Vector vecBetPointsNormalzied(Point p1, Point p2){
        Vector v;
        v.x = (p2.x-p1.x);
        v.y = (p2.y-p2.y);
        return normalizeVec(v);
    }

    Vector sumOfVectors(Vector v1, Vector v2){
        Vector v3;
        v3.x = v1.x+v2.x;
        v3.y = v1.y+v2.y;
        return normalizeVec(v3);
    }

    Vector divVecor(Vector v, int d){
        v.x /= d;
        v.y /= d;
        return normalizeVec(v);
    }

    double dMax, dMin;
    bool done;
    std::vector<AttractionPoint> attractionPoints;
    std::vector<Node> nodes;
    SpaceColonization(double maxD, double minD):dMax(maxD), dMin(minD),done(false){

    }

    void Colonize(){
        if(done) return;
        if(attractionPoints.empty() ) return;
        //for(int iAP; iAP < attractionPoints.size(); ++iAP){
        for(AttractionPoint& p : attractionPoints){
            bool discardPoint = false;
            for(int inode= 0; inode < nodes.size(); ++inode){
                Node n = nodes[inode];
                double dsq = distanceSquare(n,p);
                if( dsq < dMin*dMin){
                    attractionPoints.erase(std::remove_if(attractionPoints.begin(), attractionPoints.end(),
                                                          [&](Point t){return t.x == p.x && t.y == p.y;}), attractionPoints.end());
                    discardPoint = true;
                }else if(dsq < dMax * dMax) {
                    if(p.closestNode == -1){
                        p.closestNode = inode;
                    }else if ( dsq < distanceSquare(p, nodes[p.closestNode])){
                        p.closestNode = inode;
                    }
                }
            }

        }

        for(AttractionPoint p : attractionPoints){
            if(p.closestNode != -1){
                nodes[p.closestNode].attractors.push_back(p);
            }
        }


        std::vector<Node> newnodes;
        for(Node& n : nodes){
            if(n.attractors.size() > 0){
                Vector avg;
                for(Point p : n.attractors){
                    Vector t;
                    t.x = p.x - n.x;
                    t.y = p.y - n.y;
                    t = normalizeVec(t);
                    avg.x +=t.x;
                    avg.y +=t.y;
                }
                double angle = std::atan2(avg.y,avg.x);
                double pi = 3.141;
                Node n2;
                n2.x = n.x;
                n2.y = n.y;
                n2.direction.x = avg.x;
                n2.direction.y = avg.y;
                if(angle >= -pi && angle < -3*pi/4)
                    n2.x--;
                if(angle >= -3*pi/4 && angle < -pi/4 )
                    n2.y--;
                if(angle >= -pi/4 && angle < pi/4 )
                    n2.x++;
                else if (angle >= pi/4 && angle < 3*pi/4 )
                    n2.y++;
                else if (angle >= 3*pi/4 )
                    n2.x--;

                newnodes.push_back(n2);
            }
            n.attractors.clear();

        }
        if(newnodes.empty())
            done=true;
        for(Node& n : newnodes){
            nodes.push_back(n);
        }

    }
};

void SpaceColonizaton(){
    SpaceColonization sp(12,1);

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<> yesno(0,4);

    for(int x = 0; x < dimen; ++x){
        for(int y = 0; y<dimen; ++y){
            if(marker[x][y] != 0)
                sp.AddAttractionPoint(x,y);

        }
    }

    std::vector<Point> startPoints;
    for(int x = 0; x < dimen; ++x){
        for(int y = 0; y<dimen; ++y){
            if(map[x][y] > 0.43 && map[x][y] < 0.45){
                int r = yesno(mt);
                if(!r)
                    startPoints.push_back(Point(x,y));
            }
        }
    }

    for(auto& rm : startPoints)
        sp.AddInitalNode(rm.x, rm.y);
    for(int i = 0; i < 200; ++i)
        sp.Colonize();

    for(auto& n : sp.GetNodes()){
        marker[n.x][n.y] = 1;
    }
}


int main(){
  map.resize(dimen+1, std::vector<short>(dimen+1, 0));
  marker.resize(dimen+1, std::vector<short>(dimen+1, 0));
  DoLSystem();
}

