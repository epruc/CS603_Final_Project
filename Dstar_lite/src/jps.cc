// Copyright (c) 2018 joydeepb@cs.umass.edu
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.


#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <limits>

#include "gflags/gflags.h"
#include "eigen3/Eigen/Dense"
#include "jps.h"
#include "simple_queue.h"
#include "shared/util/timer.h"

using Eigen::Vector2i;
using std::make_pair;
using std::unordered_map;
using std::unordered_set;
using std::vector;
using std::cout;
using std::endl;

DEFINE_bool(nojps, false, "Disable JPS");
DEFINE_bool(nogui, false, "Disable Visualization");

namespace {
static const uint8_t kBlue[] = {0, 0, 255};
static const uint8_t kRed[] = {255, 0, 0};
static const uint8_t kGreen[] = {0, 255, 0};
static const uint8_t kYellow[] = {255, 255, 0};
static const uint8_t kGrey[] = {128, 128, 128};

}  // namespace

namespace astarplanner {

float Dist(const Node& v1, const Node& v2) {
  return ((v1 - v2).cast<float>().norm());
}

// Returns the path to goal, in reverse order.
void AStarPlanner::GetPath(const NodeMap& parent_map,
                           const Node& goal,
                           Path* path_ptr) {
  Node current = goal;  
  while(parent_map.find(current)!=parent_map.end()){
    path_ptr->push_back(current);
    current = parent_map.at(current);
  }
  path_ptr->push_back(current); 
}

void AStarPlanner::InitVisualization(const Map& map) {
  display_ = new cimg_library::CImgDisplay(viz_image_, "A* with JPS");
  viz_image_ = cimg_library::CImg<uint8_t>(
      map.width(), map.height(), 1, 3, 0);
  // Draw map
  for (int y = 0; y < map.height(); ++y) {
    for (int x = 0; x < map.width(); ++x) {
      if (map.Occupied(x, y)) {
        viz_image_(x, y, 0, 0) = 255;
        viz_image_(x, y, 0, 1) = 255;
        viz_image_(x, y, 0, 2) = 255;
      }
    }
  }
  acc_viz_image_ = viz_image_;
}

void AStarPlanner::Visualize(const Map& map,
                             const Node& start,
                             const Node& goal,
                             const Node& current,
                             const NodeMap& parent_map) {
  static double t_last_visualized = 0;
  static const double kMinVizInterval = 0.02;
  if (t_last_visualized > GetMonotonicTime() - kMinVizInterval) return;

  for (const auto& p : parent_map) {
    acc_viz_image_.draw_line(
      p.first.x(), p.first.y(), p.second.x(), p.second.y(), kBlue);
  }
  acc_viz_image_.draw_point(current.x(), current.y(), kYellow);

  viz_image_ = acc_viz_image_;

  viz_image_.draw_circle(current.x(), current.y(), 2, kYellow);
  viz_image_.draw_circle(start.x(), start.y(), 2, kRed);
  viz_image_.draw_circle(goal.x(), goal.y(), 2, kGreen);

  display_->display(viz_image_);
  // if (display_->is_key()) exit(0);
  t_last_visualized = GetMonotonicTime();
}

void AStarPlanner::DrawPath(const Path& path) {
  for (size_t i = 0; i + 1 < path.size(); ++i) {
    viz_image_.draw_line(
        path[i].x(),
        path[i].y(),
        path[i + 1].x(),
        path[i + 1].y(),
        kGreen);
  }
  viz_image_.draw_circle(
      path[0].x(),
      path[0].y(),
      3,
      kGreen);
  viz_image_.draw_circle(
      path.back().x(),
      path.back().y(),
      3,
      kRed);
}
std::vector<Node> gen_neigh(const Node v1, const Map& map){
  std::vector<Node> neighbors; 
  Node nextn; 
  for (int i = -1; i < 2; i=i+1){
    for(int j = -1; j < 2; j=j+1){
      if( i == 0 && j == 0){
      }else{ 
        nextn[0] = v1[0] + i; 
        nextn[1] = v1[1] + j;
        if (map.ValidNode(nextn)){
          if(!map.Occupied(nextn)){
            neighbors.push_back(nextn);
          }
        }
      }
    }
  }
  return neighbors; 
}

bool diagonal(const Node v1, const Node v2){
  int dx = v2[0]-v1[0]; 
  int dy = v2[1]-v1[1];
  if (dx == 0){
    return false; 
  }else if (dy == 0){
    return false; 
  } else { 
    return true; 
  }
}

// write a function to calculate c(x,y)
float CalcC(const Node& v1, const Node& v2,  const Map& map){
  float c = 1; 
  if(!map.ValidNode(v2)){
    c = 100000; 
  }else if(map.Occupied(v2)){
    c = 100000; 
  }else if (diagonal(v1,v2)){
    c =1.4; 
  }
  return c; 
}

float Calcg(const Node& s,
            const NodeMap& parent_map,
            const std::unordered_map<Node, float, NodeHash>& g_values,
            const Node& goal){ 
  float g = Dist(s, goal);   
  if (parent_map.find(s) != parent_map.end()) {           
    Node parent = parent_map.at(s);
    g = g_values.at(parent); 
    g = g + Dist(s,parent); 
  } 
  return g;           
}

float rhsCalc(const Node& s,
          const Node& goal, 
          const Map& map,
          const std::unordered_map<Node, float, NodeHash>& g_values){
  Node current = s; 
  std::vector<Node> neighbors = gen_neigh(s,map); 
  float rhs = 100000000; 
   for (int k = 0; k < neighbors.size(); k = k+1){ 
    Node current = neighbors[k]; 
    if (current != goal){
      float possible = Dist(s,current)+g_values.at(current); 
      rhs = std::min(rhs,possible);  
    }
  }
  return rhs; 
}

Eigen::Vector2d CalculateKeys(const Node& s,
                              const float km,
                              const Node& start, 
                              const Node& goal,
                              const NodeMap& parent_map,
                              const std::unordered_map<Node, float, NodeHash>& g_values, 
                              const Map& map){
  Eigen::Vector2d key; 
  float g = Calcg(s,parent_map,g_values,goal);
  float rhs = rhsCalc(s,goal,map,g_values); 
  float h = Dist(s,start); 
  key[0] = std::min(g, rhs + h + km); 
  key[1] = std::min(g, rhs); 
  return key; 
}
void UpdateVertex(const Node& s,
                              const float km,
                              const Node& start, 
                              const Node& goal,
                              const NodeMap& parent_map,
                              float g,
                              float rhs,
                              SimpleQueue<Node, AStarPriority>& queue,
                              const Map& map,
                              const std::unordered_map<Node, float, NodeHash>& g_values){
  if (g != rhs){
    if (g_values.find(s) == g_values.end()){
      Eigen::Vector2d knew = CalculateKeys(s,km,start, goal, parent_map, g_values,map);
      queue.Remove(s); 
      queue.Push(s,AStarPriority(knew[0],knew[1])); 
    }else{
      Eigen::Vector2d knew = CalculateKeys(s,km,start, goal, parent_map, g_values,map);
      queue.Push(s,AStarPriority(knew[0],knew[1]));  
    }
  }else if(g == rhs){
    if (g_values.find(s) != g_values.end()){
      queue.Remove(s); 
    }
  }
}

Eigen::Vector2d TopKey(SimpleQueue<Node, AStarPriority>& queue,
            std::unordered_map<Node, Eigen::Vector2d, NodeHash>& k_values){
  Node top = queue.peakTop(); 
  return k_values.at(top); 
}
bool CompareKey(std::unordered_map<Node, Eigen::Vector2d, NodeHash>& k_values,
                          SimpleQueue<Node, AStarPriority>& queue,const Node& start, 
                          const Node& goal,
                          const NodeMap& parent_map,
                          std::unordered_map<Node, float, NodeHash>& g_values,
                          std::unordered_map<Node, float, NodeHash>& rhs_values,
                          const Map& map,
                          float km,
                          const Node& s_start){
  Eigen::Vector2d topk = TopKey(queue,k_values); 
  Eigen::Vector2d startk = CalculateKeys(s_start,km,start,goal,parent_map,g_values,map); 
  if (topk[0]<startk[0]){
    if(topk[1]<startk[1]){
      return true;
    }else{
      return false;
    }
  }else{
    return false;
  }
                          }

void ComputeShortestPath(std::unordered_map<Node, Eigen::Vector2d, NodeHash>& k_values,
                          SimpleQueue<Node, AStarPriority>& queue,
                          float km,
                          const Node& start, 
                          const Node& goal,
                          const NodeMap& parent_map,
                          std::unordered_map<Node, float, NodeHash>& g_values,
                          std::unordered_map<Node, float, NodeHash>& rhs_values,
                          const Map& map,
                          const Node& s_start){
  clock_t t;
	t = clock();
  while(CompareKey(k_values,queue,start,goal,parent_map,g_values,rhs_values,map,km,s_start)){
    //cout << "h" <<endl; 
    if (rhs_values.at(s_start) > g_values.at(s_start)){ 
      break;
    }
    Node u = queue.peakTop();
    //cout << u << endl;
    //g_values[u] = Calcg(u,parent_map,g_values, goal);
    Eigen::Vector2d kold = TopKey(queue,k_values);
    Eigen::Vector2d knew = CalculateKeys(u,km,start, goal, parent_map, g_values,map);
    if((kold[1] < knew[1]) and (kold[0] < knew[0])){
      //cout << "A" << endl; 
      k_values[u] = knew;
      queue.Pop(); 
      queue.Push(u,AStarPriority(knew[0],knew[1]));
    }else if(g_values.at(u)>rhs_values.at(u)){
      //cout << "B" << endl; 
      g_values[u] = rhs_values.at(u);
      queue.Pop(); 
      rhs_values[u]= rhsCalc(u,start,map,g_values);
      UpdateVertex(u,km,start,goal,parent_map,g_values.at(u),rhs_values.at(u),queue,map,g_values);
    }else{
      //cout << "C" << endl; 
      float g_old = g_values.at(u); 
      g_values[u] = 100000000;
      std::vector<Node> neighbors = gen_neigh(u,map); 
      for (int k = 0; k < neighbors.size(); k = k+1){ 
        Node s = neighbors[k]; 
        //cout << rhs_values.at(s) << endl; 
        if (rhs_values.at(s) == CalcC(s,u,map)+g_old){
          //cout << "a" << endl; 
          if(u != goal){
            //cout << "b" << endl;
            rhs_values[u]= rhsCalc(u,start,map,g_values);
            //cout << "c" << endl; 
          }
          //cout << "e" << endl; 
        
          //cout << "f" << endl; 
          UpdateVertex(u,km,start,goal,parent_map,g_values.at(u),rhs_values.at(u),queue,map,g_values);
          //cout <<"g"<< endl; 
        }
      }
    }
    if(queue.Empty()){  
      break; 
    }
    //cout << g_values.at(u)<<endl; 
    //cout <<rhs_values.at(u)<<endl; 
    //cout <<k_values.at(u)<<endl;
  }            
}

//Dstar lite 
bool AStarPlanner::Plan(const Map& map,
                        const Node& start,
                        const Node& goal,
                        Path* path,
                        const Map& map_true) {
  int va = 0; 
  int ve = 0;
  int hp = 0;  
  clock_t t;
	t = clock();
  if (!map.ValidNode(start)) {
    printf("Invalid start location.\n");
    return false;
  }
  if (!map.ValidNode(goal)) {
    printf("Invalid goal location.\n");
    return false;
  }
  if (map.Occupied(start) || map.Occupied(goal)) {
    printf("No path possible. Start unreachable:%d goal unreachable:%d\n",
           1 - static_cast<int>(map.Occupied(start)),
           1 - static_cast<int>(map.Occupied(goal)));
    return false;
  }
  if (!FLAGS_nogui) InitVisualization(map);
  Node s_last = start; 
  Node s_start = start;  
  // Initialize parent map.
  parent_map_.clear();
  // Clear all G values.
  g_values_.clear();
  // Clear all rhs values.
  rhs_values_.clear();
  // Clear all rhs values.
  k_values_.clear();
  // Clear all k values.
  cost_true_.clear(); 
  cost_.clear();
  for(int i = 0 ; i < map.width(); i = i + 1){
    for(int j = 0; j< map.height(); j = j+1){
      Node current; 
      current[0] = i;
      current[1] = j;
      rhs_values_[current] = 100000000;
      g_values_[current] = 100000000; 
    }
  } 
  // Initialize an empty priority queue.
  SimpleQueue<Node, AStarPriority> queue;
  // Add start to priority queue.
  const float goal_heuristic = Dist(start,goal); // TODO
  queue.Push(goal, AStarPriority(goal_heuristic, 0));
  va = va +1; 
  // Clear the closed set.
  closed_set_.clear();
  //Initilization 
  float km = 0;
  rhs_values_[goal] = 0;
  g_values_[goal] = 0;  
  k_values_[goal] = CalculateKeys(goal,km,start, goal, parent_map_, g_values_,map);
  ComputeShortestPath(k_values_,queue,km,start, goal, parent_map_, g_values_, rhs_values_,map,s_start);
  hp = hp + 1; 
  //cout << "here" << endl; 
  //Visualize(map, start, goal,s_start, parent_map_)
  while(s_start != goal){
    //cout << s_start << endl; 
    path->push_back(s_start); 
    // if (rhs_values_[s_start] ==100000000){
    //   cout <<"No possible Path"<< endl; 
    // }
    std::vector<Node> neighbors = gen_neigh(s_start, map); 
    va = va +1; 
    ve = ve +1;
    for (int k = 0; k < neighbors.size(); k = k+1){ 
      Node current = neighbors[k]; 
      va = va +1; 
      cost_[current] = CalcC(s_start,current,map);
      cost_true_[current] = CalcC(s_start,current,map_true);
    } 
    float min = 100000000;
    Node next = s_start; 
    va = va +1; 
    ve = ve +1;
    for (int l = 0; l < neighbors.size(); l = l+1){ 
      Node pos = neighbors[l]; 
      va = va +1; 
      float cost = cost_true_[pos] + Calcg(pos,parent_map_,g_values_,goal); 
      if (cost < min){
        next = pos;
        min = cost; 
      }
    }
    //cout << "here" << endl; 
    s_start = next; 
    std::vector<Node> neighbor = gen_neigh(s_start, map); 
    float flag = 0;
    va = va +1; 
    ve = ve +1;
    for (int m = 0; m < neighbor.size(); m = m+1){ 
      Node current = neighbor[m]; 
      va = va +1; 
      cost_[current] = CalcC(s_start,current,map);
      cost_true_[current] = CalcC(s_start,current,map_true);
      if (cost_.at(current) != cost_true_.at(current)){
        flag = 1; 
        if (cost_.at(current) > cost_true_.at(current)){
          if(s_start != goal){
            rhs_values_[s_start] = std::min(rhsCalc(s_start,goal,map_true,g_values_),CalcC(s_start,current,map_true)+Calcg(current,parent_map_,g_values_,goal));
          }
        }else if (rhs_values_[s_start] == CalcC(s_start,current,map) + Calcg(current,parent_map_,g_values_,goal)){
          if(s_start != goal){
            rhs_values_[s_start] = rhsCalc(s_start,goal,map_true,g_values_);
          }
        }
        UpdateVertex(s_start,km,start,goal,parent_map_,g_values_.at(start),rhs_values_.at(start),queue,map,g_values_);
        hp = hp + 1;
      }
    } 
    
    if (flag == 1){
      km = km + Dist(s_last,goal); 
      s_last = s_start; 
      if (!queue.Empty()){
        ComputeShortestPath(k_values_,queue,km,start, goal, parent_map_, g_values_, rhs_values_,map,s_start);
        hp = hp +1;
      }
      //Visualize(map, start, goal,s_start, parent_map_);
    }
  }
  t = clock() - t;
	cout << "time: " << t*1.0/CLOCKS_PER_SEC << " seconds" << endl;
  cout <<"ve: " << ve << endl;
  cout << "va: " << va << endl;
  cout << "hp: " << hp << endl; 
  // To visualize the path found:
  if (!FLAGS_nogui) {
    DrawPath(*path);
    display_->display(viz_image_);
    while (!display_->is_closed() && !display_->is_key()) display_->wait();
  }

  return false;
}

}  // namespace astarplanner

