//
// Created by lchavarr on 4/17/16.
//


#include "solctra_sequential.h"
#include <omp.h>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstring>
#include <string>
#include <sstream>

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


void RK4(const GlobalData& data, const std::string& output, const cartesian& start_point, const int steps, const double& step_size, const int particle, const int mode, const int print_type)
{
    Coil rmi[TOTAL_OF_COILS];
    Coil rmf[TOTAL_OF_COILS];
    initializeGlobals(rmi, rmf);
    cartesian p0;
    cartesian p1 = {0, 0, 0};
    cartesian p2 = {0, 0, 0};
    cartesian p3 = {0, 0, 0};
    cartesian K1;
    cartesian K2;
    cartesian K3;
    cartesian K4;
    cartesian K5;
    cartesian Ovect = {0, 0, 0};
    cartesian p = {0, 0, 0};
    cartesian r_vector;
    double norm_temp;
    double r_radius;

    /*FILE* handler;
    std::string file_name = output +  "/path" + getZeroPadded(particle) + ".txt";
    handler = fopen(file_name.c_str(), "w");
    if(nullptr == handler)
    {
        printf("Unable to open file=[%s]. Nothing to do\n", file_name.c_str());
        exit(0);
    }

    if (print_type == 0)
    {
	fprintf(handler, "%e\t%e\t%e\n", start_point.x, start_point.y, start_point.z);
    }

    if (print_type == 1)
    {
	fprintf(handler, "x,y,z\n");
	fprintf(handler, "%e,%e,%e\n", start_point.x, start_point.y, start_point.z);

    }*/


    p0 = start_point;
    //const double steps_inverse = static_cast<double>(1) / steps;
    //const int onePercent = static_cast<int>(steps / 100);
    const double half = 1.0 / 2.0;

    for (int i = 1; i < steps; i++)
    {
        K1 = magnetic_field(rmi, rmf, data, p0);
        //printf("After magnetic fields.\n");
        norm_temp = 1.0 / norm_of(K1);
        K1.x = ( K1.x * norm_temp ) * step_size;
        K1.y = ( K1.y * norm_temp ) * step_size;
        K1.z = ( K1.z * norm_temp ) * step_size;
        p1.x = ( K1.x * half ) + p0.x;
        p1.y = ( K1.y * half ) + p0.y;
        p1.z = ( K1.z * half ) + p0.z;

        K2 = magnetic_field(rmi, rmf, data, p1);
        norm_temp = 1.0 / norm_of(K2);
        K2.x = ( K2.x * norm_temp ) * step_size;
        K2.y = ( K2.y * norm_temp ) * step_size;
        K2.z = ( K2.z * norm_temp ) * step_size;
        p2.x = ( K2.x * half ) + p0.x;
        p2.y = ( K2.y * half ) + p0.y;
        p2.z = ( K2.z * half ) + p0.z;

        K3 = magnetic_field(rmi, rmf, data, p2);
        norm_temp = 1.0 / norm_of(K3);
        K3.x = ( K3.x * norm_temp ) * step_size;
        K3.y = ( K3.y * norm_temp ) * step_size;
        K3.z = ( K3.z * norm_temp ) * step_size;
        p3.x = K3.x + p0.x;
        p3.y = K3.y + p0.y;
        p3.z = K3.z + p0.z;

        K4 = magnetic_field(rmi, rmf, data, p3);
        norm_temp = 1.0 / norm_of(K4);
        K4.x = ( K4.x * norm_temp ) * step_size;
        K4.y = ( K4.y * norm_temp ) * step_size;
        K4.z = ( K4.z * norm_temp ) * step_size;
        p0.x = p0.x + (( K1.x + 2 * K2.x + 2 * K3.x + K4.x ) / 6 );
        p0.y = p0.y + (( K1.y + 2 * K2.y + 2 * K3.y + K4.y ) / 6 );
        p0.z = p0.z + (( K1.z + 2 * K2.z + 2 * K3.z + K4.z ) / 6 );


	//K5 = magnetic_field(rmi, rmf, data, p0);


	/*if (print_type == 0)
	{
		//fprintf(handler, "%e\t%e\t%e\t%e\n", p0.x, p0.y, p0.z,norm_of(K5));
		 fprintf(handler, "%e\t%e\t%e\n", p0.x, p0.y, p0.z);
	}

	if (print_type == 1)
	{
		//fprintf(handler, "%e,%e,%e,%e\n", p0.x, p0.y, p0.z, norm_of(K5));
		 fprintf(handler, "%e,%e,%e\n", p0.x, p0.y, p0.z);
	}*/

        if (mode == 1)
        {
            p.x = p0.x;
            p.y = p0.y;
            Ovect.x = ( p.x / norm_of(p)) * 0.2381; //// Origen vector
            Ovect.y = ( p.y / norm_of(p)) * 0.2381;
            Ovect.z = 0;
            r_vector.x = p0.x - Ovect.x;
            r_vector.y = p0.y - Ovect.y;
            r_vector.z = p0.z - Ovect.z;
            r_radius = norm_of(r_vector);
            if (r_radius > 0.0944165)
            {
                //fprintf(handler, "%e,%e,%e\n", r_radius, 0.0, 0.0);
                break;
            }
        }
    }
    //fclose(handler);
    finishGlobal(rmi,rmf);
}


void getMagneticProfile(const GlobalData& data, const int num_points, const int phi_angle, const std::string& output, const int dimension){

 	//Prepare parameters for magnetic_field function: rmi, rmf
	Coil rmi[TOTAL_OF_COILS];
	Coil rmf[TOTAL_OF_COILS];
	Coil observation_points;
	cartesian point={0,0,0};
	cartesian B_point;
	const double major_R = 0.2381;
	const double minor_r = 0.0944165;
	double width;
 	double radians = phi_angle*PI/180.0;

    initializeGlobals(rmi, rmf);

    //Prepare output file
    FILE* handler;
    std::string file_name = output + "/magnetic_field.txt";
    std::cout << "Before Handler open" << std::endl;
    handler = fopen(file_name.c_str(), "w");
    std::cout << "After handler open" << std::endl;
    if(nullptr == handler)
    {
        printf("Unable to open file=[%s]. Nothing to do\n", file_name.c_str());
        exit(0);
    }
   
    fprintf(handler, "x,y,z,|B|\n");
   	
    if(dimension == 1){
        width = (2*minor_r)/num_points;
        observation_points.x = static_cast<double*>(malloc(sizeof(double) * num_points));
        observation_points.y = static_cast<double*>(malloc(sizeof(double) * num_points));
        observation_points.z = static_cast<double*>(malloc(sizeof(double) * num_points));
        //Generate observation points at phi_angle plane
        for(int i=0; i<num_points; i++){
            observation_points.x[i] = ((major_R-minor_r+(width*i))+minor_r*cos(PI/2))*cos(radians);
            observation_points.y[i] = ((major_R-minor_r+(width*i))+minor_r*cos(PI/2))*sin(radians);
            observation_points.z[i] = 0.0;
        }
        for(int i=0;i<num_points;i++)
        {
            point.x = observation_points.x[i];
            point.y = observation_points.y[i];
            point.z = observation_points.z[i];
            B_point = magnetic_field(rmi,rmf,data,point);
            fprintf(handler, "%e,%e,%e,%e\n", point.x, point.y, point.z,norm_of(B_point));
        }
    }else if(dimension == 2){
        std::cout << "Entrando a dimension 2" << std::endl;
	width = minor_r/num_points;
        observation_points.x = static_cast<double*>(malloc(sizeof(double) * (num_points*360)));
        observation_points.y = static_cast<double*>(malloc(sizeof(double) * (num_points*360)));
        observation_points.z = static_cast<double*>(malloc(sizeof(double) * (num_points*360)));
        std::cout << "Inicializando puntos de observacion" << std::endl;
	for(int i=0; i<360; i++){
            for(int j=0; j<num_points; j++){
		    observation_points.x[((num_points*i)+j)] = (major_R+((width*j)*sin(i*(PI/180))))*cos(radians);
        	    observation_points.y[((num_points*i)+j)] = ((width*j)*cos(i*PI/180));
           	    observation_points.z[((num_points*i)+j)] = (major_R+(width*j)* sin(i*PI/180))*sin(radians);
            }
        }
	std::cout << "Inicializacion finalizada" << std::endl;

        for(int i=0;i<num_points*360;i++)
        {
            point.x = observation_points.x[i];
            point.y = observation_points.y[i];
            point.z = observation_points.z[i];
            B_point = magnetic_field(rmi,rmf,data,point);
            fprintf(handler, "%e,%e,%e,%e\n", point.x, point.y, point.z,norm_of(B_point));
        }
	std::cout << "Campo calculado" << std::endl;

    }
     //For each observation point call magnetic_field
	fclose(handler);
	free(observation_points.x);
	free(observation_points.y);
	free(observation_points.z);
	finishGlobal(rmi, rmf);
}


void runParticles(const GlobalData& data, const std::string& output, const Coil& particles, const int length, const int steps, const double& step_size, const int mode, const int print_type)
{
    cartesian A={0,0,0};
    for(int i=0; i < length ; ++i)
    {
        A.x = particles.x[i];
        A.y = particles.y[i];
        A.z = particles.z[i];
        A.print();
        RK4(data, output, A,steps,step_size,i, mode, print_type);

    }
}
