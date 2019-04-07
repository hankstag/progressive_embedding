

#include <queue>
#include <algorithm>
#include <fstream>
#include <iostream>

#include "is_simple_polygon.h"
#include <igl/copyleft/cgal/orient2D.h>
#include <igl/copyleft/cgal/segment_segment_intersect.h>


typedef std::tuple<double,double,int> Point;
typedef std::tuple<Point,Point,int> Segment;
typedef std::tuple<Point,int,bool> Event; // (point, segment_id, enter/leave)

bool is_simple_result = true;

struct Order{
    bool operator() (const Segment& a, const Segment& b) const {
        int id1 = std::get<2>(a);
        int id2 = std::get<2>(b);
        if(id1 == id2) return false;
        // TODO: proper comparison
        // notation: [a0 a1] [b0 b1]
        // ordering rule: 
        // | if a0.x != b0.x : bool r = left_is_higher(seg1, seg2)
        // | if a0.x == b0.x :                                    
        // |      if a0.y == b0.y:                                 
        // |          switch orient(a0,b0,b1):                     
        // |          case 1:  k   
        // |          case -1: k
        // |          case 0: not simple
        // |      else 
        // |          return a0.y < b0.y
        auto left_is_higher = [&](const Segment& seg1, const Segment& seg2){
            double al[2] = {std::get<0>(std::get<0>(seg1)), std::get<1>(std::get<0>(seg1))};
            double ar[2] = {std::get<0>(std::get<1>(seg1)), std::get<1>(std::get<1>(seg1))};
            double bl[2] = {std::get<0>(std::get<0>(seg2)), std::get<1>(std::get<0>(seg2))};
            double br[2] = {std::get<0>(std::get<1>(seg2)), std::get<1>(std::get<1>(seg2))};
            switch(igl::copyleft::cgal::orient2D(al,ar,bl)){
                case -1: return false;break;
                case 1 : return true ;break;
                case 0 : if(std::get<2>(std::get<0>(seg2)) != std::get<2>(std::get<0>(seg1)) &&
                            std::get<2>(std::get<0>(seg2)) != std::get<2>(std::get<1>(seg1)))
                            // b1 and a0 a1 are not the same point
                            is_simple_result = false;
                         else return (al[1]<bl[1]);
            }
            return false;
        };
        Point a0 = std::get<0>(a);
        Point b0 = std::get<0>(b);
        if(std::get<0>(a0) != std::get<0>(b0)){
            return (std::get<0>(a0) > std::get<0>(b0)) ? !left_is_higher(b,a) : left_is_higher(a,b);
        }else{
            if(std::get<1>(a0) == std::get<1>(b0)){ // or use vertex id
                double al[2] = {std::get<0>(std::get<0>(a)), std::get<1>(std::get<0>(a))};
                double ar[2] = {std::get<0>(std::get<1>(a)), std::get<1>(std::get<1>(a))};
                double bl[2] = {std::get<0>(std::get<0>(b)), std::get<1>(std::get<0>(b))};
                double br[2] = {std::get<0>(std::get<1>(b)), std::get<1>(std::get<1>(b))};
                switch(igl::copyleft::cgal::orient2D(al,ar,br)){
                    case -1: return false;break;
                    case 1 : return true ;break;
                    case 0 : is_simple_result = false;
                }
                return false;
            }else 
                return std::get<1>(a0) < std::get<1>(b0);
        }

    }
};
using SweepList = std::set<Segment,Order>;

// check whether Segment a and Segment b intersect each other
bool disjoint_segment_intersect(
    Segment s1,
    Segment s2,
    int n
){
    int id1 = std::get<2>(s1);
    int id2 = std::get<2>(s2);
    if(std::abs(id1-id2)==1 || std::abs(id1-id2)==n-1){
        // if s1 and s2 are exactly the same
        if(std::get<0>(s1) == std::get<0>(s2) &&
           std::get<1>(s1) == std::get<1>(s2))
            return true;
        else
            return false;
    }
    double a[2] = {std::get<0>(std::get<0>(s1)),std::get<1>(std::get<0>(s1))};
    double b[2] = {std::get<0>(std::get<1>(s1)),std::get<1>(std::get<1>(s1))};
    double c[2] = {std::get<0>(std::get<0>(s2)),std::get<1>(std::get<0>(s2))};
    double d[2] = {std::get<0>(std::get<1>(s2)),std::get<1>(std::get<1>(s2))};
    return igl::copyleft::cgal::segment_segment_intersect(a,b,c,d,0.0);
}

bool is_simple_polygon(const Eigen::MatrixXd& P){
    // P: [a, b, c, d,...]
    // change it to the form of edges
    // E: [a, b; b, c; c, d;...]
    is_simple_result = true;
    // build segments list
    std::vector<Segment> segments;
    for(int i=0;i<P.rows();i++){
        int i_1 = (i+1)%P.rows();
        Point p1 = Point(P(i,  0),P(i  ,1),i  );
        Point p2 = Point(P(i_1,0),P(i_1,1),i_1);
        if(p1>p2)
            segments.push_back(Segment(p2,p1,segments.size()));
        else 
            segments.push_back(Segment(p1,p2,segments.size()));
    }

    // build event list (sort segments)
    auto later = [](const Event& a, const Event& b){
        if(std::get<0>(a) == std::get<0>(b)){
            if(std::get<2>(a) == true && std::get<2>(b) == false)
                return false;
            else if(std::get<2>(a) == false && std::get<2>(b) == true)
                return true;
            else 
                return a > b;
        }else
            return a > b;
    };
    std::priority_queue<Event,std::vector<Event>,decltype(later)> Q(later);
    for(int i=0;i<segments.size();i++){
        Q.push(Event(std::get<0>(segments[i]),i,true ));
        Q.push(Event(std::get<1>(segments[i]),i,false));
    }
    SweepList sl;
    int n = segments.size();
    // start line sweeping
    while(!Q.empty()){
        Event evt = Q.top();
        Q.pop();
        Segment seg=segments[std::get<1>(evt)];
        bool is_enter = std::get<2>(evt);
        int seg_id = std::get<2>(seg);
        if(is_enter){
            sl.insert(seg);
            auto pos = sl.find(seg);
            auto prev = (pos == sl.begin()) ? pos: std::prev(pos);
            auto next = std::next(pos);
            if((pos  != sl.begin() && disjoint_segment_intersect(*pos,*prev,n)) ||
               (next != sl.end()   && disjoint_segment_intersect(*pos,*next,n)))
                return false;
        }else{
            auto pos = sl.find(seg);
            auto prev = (pos == sl.begin()) ? pos : std::prev(pos);
            auto next = (pos == sl.end()) ? pos : std::next(pos);
            if(pos != sl.begin() && next != sl.end() && disjoint_segment_intersect(*prev,*next,n))
                return false;
            sl.erase(seg);
        }
    }
    return is_simple_result;
}

