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

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "gflags/gflags.h"

#include "jps.h"
#include "util/random.h"
#include "util/timer.h"

using astarplanner::Map;
using astarplanner::Node;
using astarplanner::Path;
using astarplanner::AStarPlanner;

void Plan(const Map& map, const Node& start, const Node& goal,const Map& map_true) {
  AStarPlanner planner(map.width());
  Path path;
  planner.Plan(map, start, goal, &path, map_true);
}

int main(int num_args, char* args[]) {
  // google::InitGoogleLogging(args[0]);
  google::ParseCommandLineFlags(&num_args, &args, true);
  if (num_args < 2) {
    printf("Usage: ./bin/jps map.png [goal_x goal y [start x start y]]\n");
    return 1;
  }
  printf("Planning on map image '%s'\n", args[1]);
  Map map(args[1]);
  Map map_true(args[6]);
  Node start(15, 15);
  Node goal(map.width() - 10, map.height() / 2);
  if (num_args > 3) {
    goal.x() = atoi(args[2]);
    goal.y() = atoi(args[3]);
  }
  if (num_args > 5) {
    start.x() = atoi(args[4]);
    start.y() = atoi(args[5]);
  }
  printf("Start: %d,%d Goal: %d,%d\n",
          start.x(), start.y(),
          goal.x(), goal.y());
  Plan(map, start, goal, map_true);
  return 0;
}
