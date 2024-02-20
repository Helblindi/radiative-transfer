#include "param.h"

namespace parameter 
{
   template <>
   bool parameter::get<bool>(std::string key, bool value) {
   if (!contains(key)) {
      return value;
   }
   if ("yes" == params[key] || "Yes" == params[key] ||
       "true" == params[key] || "True" == params[key]) {
         std::cout << "Found " << key << " to be true.\n";
      return true;
   } else {
      return false;
   }
   }

   template <>
   bool parameter::get(std::string key) {
   check_key(key);
   return get<bool>(key, false);
   }

   template <>
   int parameter::get(std::string key, int value) {
   if (!contains(key)) {
      return value;
   }
   return std::stoi(params[key]);
   }

   template <>
   int parameter::get(std::string key) {
   check_key(key);
   return get<int>(key, 0);
   }

   template <>
   double parameter::get(std::string key, double value) {
   if (!contains(key)) {
      return value;
   }
   return std::stod(params[key]);
   }

   template <>
   double parameter::get(std::string key) {
   check_key(key);
   return get<double>(key, 0.0);
   }

   // String addition
   template <>
   std::string parameter::get(std::string key, std::string value) {
   if (!contains(key)) {
      return value;
   }
   return params[key];
   }

   template <>
   std::string parameter::get(std::string key) {
   check_key(key);
   return get<std::string>(key, "");
   }
} // end namespace param