#pragma once

#include <iostream>

#include <chrono>

class MachineTimer
{
public:
    MachineTimer(void)
    {
        elapsed = 0;

        start = std::chrono::high_resolution_clock::now();
    };

    void Start(){start = std::chrono::high_resolution_clock::now();}

    void Stop()
    {
        elapsed += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
    }

    void PrintElapsed()
    {
        std::cout << "Time taken = " << elapsed << '\n';
    }

    void ResetPrintElapsed(){Stop();PrintElapsed();Start();}

    void Reset(){elapsed = 0;}

    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    double elapsed;
};
