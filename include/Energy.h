#pragma once

struct Group
{
   double eg_min, eg_max;

   Group(double group_min, double group_max)
     : eg_min(group_min), eg_max(group_max)
   {}
};
