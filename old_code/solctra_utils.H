//
// Created by lchavarr on 4/17/16.
//

#ifndef SOLCTRA_UTILS_H
#ifndef __INTEL_COMPILER
#define SOLCTRA_UTILS_H


#include <string>
#include <cmath>
#include <sstream>
#include <cstdio>

#define     PI      3.141592654
#define     miu     1.2566e-06
#define     I       -4350
#define ALIGNMENT_SIZE 64

#ifdef KNL
#define GRADES_PER_PAGE ALIGNMENT_SIZE  * KNL / sizeof(double)
#else
#define GRADES_PER_PAGE ALIGNMENT_SIZE / sizeof(double)
#endif

#define TOTAL_OF_GRADES 360
#define TOTAL_OF_GRADES_PADDED 384
#define TOTAL_OF_COILS 12
//#define PATH_TO_RESOURCES "resources"

struct cartesian
{
    double x, y, z;
    void print()
    {
        printf("X=[%e]. Y=[%e]. Z=[%e].\n", x, y, z);
    }
};

struct Coil
{
    double* x;
    double* y;
    double* z;

};

struct GlobalData
{
    Coil coils;
    Coil e_roof;
    double* leng_segment;
};

void* _mm_malloc(size_t size, size_t alignment);

void _mm_free(void* pointer);

void loadFile(double* x, double* y, double* z, const int length, const std::string& path);

double getCurrentTime();

void createDirectoryIfNotExists(const std::string& path);

bool directoryExists(const std::string& path);

std::string getZeroPadded(const int num);

void load_coil_data(double* x, double* y, double* z, 
                    const std::string& path);

void e_roof(GlobalData& data);

void R_vectors(const Coil& coil, 
                const cartesian& point, 
                Coil* Rmi, Coil* Rmf);

cartesian magnetic_field(Coil* rmi, Coil* rmf, 
                        const GlobalData& data, 
                        const cartesian& point);

cartesian magField(Coil* rmi, Coil* rmf, 
                        const GlobalData& data, 
                        const cartesian& point);

void initializeGlobals(Coil* rmi, Coil* rmf);

void finishGlobal(Coil* rmi, Coil* rmf);

inline double norm_of(const cartesian& vec){
    return sqrt(( vec.x * vec.x ) + ( vec.y * vec.y ) + ( vec.z * vec.z ));
}

GlobalData initialize_coils();

void free_coils(GlobalData& data);

void* _mm_malloc(size_t size, size_t alignment);
void _mm_free(void* pointer);

void loadFile(double* x, double* y, double* z, 
                const int length, const std::string& path);



#endif //__INTEL_COMPILER
#endif // SOLCTRA_UTILS_H
