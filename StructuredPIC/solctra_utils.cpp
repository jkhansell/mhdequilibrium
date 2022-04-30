//
// Created by lchavarr on 4/17/16.
//

#include "solctra_utils.h"
#include <cstdio>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstring>
#include <string>
#include <sstream>


void* _mm_malloc(size_t size, size_t /*alignment*/)
{
    return malloc(size);
}
void _mm_free(void* pointer)
{
    free(pointer);
}

void initializeGlobals(Coil* rmi, Coil* rmf)
{
    for(unsigned int i = 0 ; i < TOTAL_OF_COILS ; ++i)
    {
        rmi[i].x = static_cast<double*>(_mm_malloc(sizeof(double) * (TOTAL_OF_GRADES + 1), ALIGNMENT_SIZE));
        rmi[i].y = static_cast<double*>(_mm_malloc(sizeof(double) * (TOTAL_OF_GRADES + 1), ALIGNMENT_SIZE));
        rmi[i].z = static_cast<double*>(_mm_malloc(sizeof(double) * (TOTAL_OF_GRADES + 1), ALIGNMENT_SIZE));
        rmf[i].x = static_cast<double*>(_mm_malloc(sizeof(double) * (TOTAL_OF_GRADES + 1), ALIGNMENT_SIZE));
        rmf[i].y = static_cast<double*>(_mm_malloc(sizeof(double) * (TOTAL_OF_GRADES + 1), ALIGNMENT_SIZE));
        rmf[i].z = static_cast<double*>(_mm_malloc(sizeof(double) * (TOTAL_OF_GRADES + 1), ALIGNMENT_SIZE));
    }
}

void finishGlobal(Coil* rmi, Coil* rmf)
{
    for(unsigned int i = 0 ; i < TOTAL_OF_COILS ; ++i)
    {
        _mm_free(rmi[i].x);
        _mm_free(rmi[i].y);
        _mm_free(rmi[i].z);
        _mm_free(rmf[i].x);
        _mm_free(rmf[i].y);
        _mm_free(rmf[i].z);
    }
}
void load_coil_data(double* x, double* y, double* z, const std::string& path)
{
    for (int num = 0; num < TOTAL_OF_COILS; num++)
    {
	std::ostringstream convert;
	convert << num;
	std::string value = convert.str();
	std::string tmp =  path + "/Bobina"+value+"m.txt";
        loadFile(&(x[num * TOTAL_OF_GRADES_PADDED]), &(y[num * TOTAL_OF_GRADES_PADDED]), &(z[num * TOTAL_OF_GRADES_PADDED]), TOTAL_OF_GRADES + 1, tmp);
    }
}

void e_roof(GlobalData& data)
{
    cartesian segment;
    for (int j = 0; j < TOTAL_OF_COILS; j++)
    {
        const int base = j * TOTAL_OF_GRADES_PADDED;
        for (int i = 0; i < TOTAL_OF_GRADES; i++)
        {
            segment.x = ( data.coils.x[base + i + 1] ) - ( data.coils.x[base + i] );
            segment.y = ( data.coils.y[base + i + 1] ) - ( data.coils.y[base + i] );
            segment.z = ( data.coils.z[base + i + 1] ) - ( data.coils.z[base + i] );
            data.leng_segment[base + i] = norm_of(segment);
            const double leng_segment_inverted = 1.0 / data.leng_segment[base + i];
            data.e_roof.x[base + i] = segment.x * leng_segment_inverted;
            data.e_roof.y[base + i] = segment.y * leng_segment_inverted;
            data.e_roof.z[base + i] = segment.z * leng_segment_inverted;
        }
    }
}

void R_vectors(const Coil& coil, const cartesian& point, Coil* Rmi, Coil* Rmf)
{
    for(unsigned int i = 0 ; i < TOTAL_OF_COILS ; ++i)
    {
        const int base = i * TOTAL_OF_GRADES_PADDED;
        double* x = &coil.x[base];
        double* y = &coil.y[base];
        double* z = &coil.z[base];
        for (int j = 0; j < TOTAL_OF_GRADES; j++)
        {
            Rmi[i].x[j] = point.x - x[j];
            Rmi[i].y[j] = point.y - y[j];
            Rmi[i].z[j] = point.z - z[j];
        }
        for (int j = 0; j < TOTAL_OF_GRADES; j++)
        {
            Rmf[i].x[j] = point.x - x[j + 1];
            Rmf[i].y[j] = point.y - y[j + 1];
            Rmf[i].z[j] = point.z - z[j + 1];
        }
    }
}

cartesian magnetic_field(Coil* rmi, Coil* rmf, const GlobalData& data, const cartesian& point)
{
    const double multiplier = ( miu * I ) / ( 4 * PI );

    double Bx = 0;
    double By = 0;
    double Bz = 0;

    for (int i = 0; i < TOTAL_OF_COILS; i++)
    {
        for (int jj = 0; jj < TOTAL_OF_GRADES; jj += GRADES_PER_PAGE)
        {
            const unsigned final = (TOTAL_OF_GRADES < jj + GRADES_PER_PAGE) ? TOTAL_OF_GRADES : jj + GRADES_PER_PAGE;
            const int base = i * TOTAL_OF_GRADES_PADDED;
            double* x = &data.coils.x[base];
            double* y = &data.coils.y[base];
            double* z = &data.coils.z[base];
            for (int j = jj; j < final ; ++j)
            {
                rmi[i].x[j] = point.x - x[j];
                rmi[i].y[j] = point.y - y[j];
                rmi[i].z[j] = point.z - z[j];
                rmf[i].x[j] = point.x - x[j + 1];
                rmf[i].y[j] = point.y - y[j + 1];
                rmf[i].z[j] = point.z - z[j + 1];
            }
            for (int j = jj; j < final ; ++j)
            {
               const double norm_Rmi = sqrt((( rmi[i].x[j] * rmi[i].x[j] ) + ( rmi[i].y[j] * rmi[i].y[j] ) +
                                              ( rmi[i].z[j] * rmi[i].z[j] )));
               const double norm_Rmf = sqrt((( rmf[i].x[j] * rmf[i].x[j] ) + ( rmf[i].y[j] * rmf[i].y[j] ) +
                                              ( rmf[i].z[j] * rmf[i].z[j] )));

                //firts vector of cross product in equation 8
                cartesian U;
                U.x = multiplier * data.e_roof.x[base + j];
                U.y = multiplier * data.e_roof.y[base + j];
                U.z = multiplier * data.e_roof.z[base + j];

                //second vector of cross product in equation 8
                const double C = (
                        (( 2 * ( data.leng_segment[base + j] ) * ( norm_Rmi + norm_Rmf )) /
                          ( norm_Rmi * norm_Rmf )) *
                         (( 1 ) / (( norm_Rmi + norm_Rmf ) * ( norm_Rmi + norm_Rmf ) -
                                  data.leng_segment[base + j] * data.leng_segment[base + j] )));

                cartesian V;
                V.x = rmi[i].x[j] * C;
                V.y = rmi[i].y[j] * C;
                V.z = rmi[i].z[j] * C;

                //cross product in equation 8
                Bx = Bx + (( U.y * V.z ) - ( U.z * V.y ));
                By = By - (( U.x * V.z ) - ( U.z * V.x ));
                Bz = Bz + (( U.x * V.y ) - ( U.y * V.x ));
            }
        }
    }

    cartesian B = {0.0, 0.0, 0.0};
    B.x = Bx;
    B.y = By;
    B.z = Bz;

    return B;
}

GlobalData initialize_coils(){
    /*
    Coil data initializer

    Obtains data from the resources folder which contains the coordinates of SCR-1's modular coils.
    */

    std::string resourcePath = "./resources";   //Coil directory

    const size_t sizeToAllocate = sizeof(double) * TOTAL_OF_GRADES_PADDED * TOTAL_OF_COILS;     //size of allocation

    GlobalData data;    //coil data structure

    // memory allocation
    data.coils.x = static_cast<double*>(_mm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    data.coils.y = static_cast<double*>(_mm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    data.coils.z = static_cast<double*>(_mm_malloc(sizeToAllocate, ALIGNMENT_SIZE));

    // coil data loading
    load_coil_data(data.coils.x, data.coils.y, data.coils.z, resourcePath);
    
    // e_roof memory allocation
    data.e_roof.x = static_cast<double*>(_mm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    data.e_roof.y = static_cast<double*>(_mm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    data.e_roof.z = static_cast<double*>(_mm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    data.leng_segment = static_cast<double*>(_mm_malloc(sizeToAllocate, ALIGNMENT_SIZE));

    // e_roof calculation
    e_roof(data);

    return data;  

}

// free coil/e_roof memory
void free_coils(GlobalData& data){
    _mm_free(data.coils.x);
    _mm_free(data.coils.y);
    _mm_free(data.coils.z);
    _mm_free(data.e_roof.x);
    _mm_free(data.e_roof.y);
    _mm_free(data.e_roof.z);
    _mm_free(data.leng_segment);
}

void loadFile(double* x, double* y, double* z, const int length, const std::string& path)
{
    FILE* file_buff;
    //Open file
    file_buff = fopen(path.c_str(), "r");
    if (file_buff == nullptr)
    {
        printf("Error al abrir archivo \n");
    }
    else
    {
        double localX, localY, localZ;
        printf("Loading %s with length=%d\n", path.c_str(), length);
        for (int point = 0; point < length; point++)
        {
            fscanf(file_buff, "%le %le %le", &localX, &localY, &localZ);
            x[point] = localX;
            y[point] = localY;
            z[point] = localZ;
        }
        fclose(file_buff);
    }
}
