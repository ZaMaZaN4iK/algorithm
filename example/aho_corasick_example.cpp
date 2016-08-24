/*
   Copyright (c) Alexander Zaitsev <zamazan4ik@gmail.by>, 2016

   Distributed under the Boost Software License, Version 1.0. (See accompanying
   file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

    For more information, see http://www.boost.org
*/

#include <vector>
#include <string>
#include <iostream>

#include <boost/algorithm/searching/aho_corasick.hpp>


int main()
{
    std::vector<std::string> pat({"228", "he", "is", "1488", "she", "his", "322", "her",
                                  "h", "hishera", "azaza"});
    std::string corp = "hisher";
    std::vector<std::pair<std::string::const_iterator, std::string::const_iterator>> out;

    bool result = boost::algorithm::aho_corasick_map<char>(corp.begin(), corp.end(), pat.begin(), pat.end(), 
					     [&out](std::string::const_iterator begin, std::string::const_iterator end) -> bool
					     { out.push_back({begin, end}); return true; });

    std::cout << result << std::endl;
    for(const auto& val: out)
    {
        auto begin = val.first;
        auto end = val.second;
        while (begin != end)
        {
            std::cout << *begin;
            ++begin;
        }
        std::cout << std::endl;
    }
    return 0;
}