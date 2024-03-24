//
// Created by andi on 30.05.23.
//

#include "Polymorphism.h"
#include <iostream>
#include <chrono>
#include <stdio.h>
#include <vector>
#include <numeric>
#include <boost/asio/ip/host_name.hpp>
#include "../Header/Helpfunctions and Classes.h"

using namespace std;


int main() {
    std::string timeString = "1:00:00:00";
    int minutes = timeStringToMinutes(timeString);
    std::cout << "Number of minutes: " << minutes << std::endl;

    string teststring = exec("echo test");
    cout << "printing teststring: " << teststring;

    cout << timeStringToMinutes("12:00") << endl;
    cout << timeStringToMinutes("1:12:00") << endl;

    cout << get_remaining_minutes() << endl;

    return 0;
}