
#include "solctra_sequential.h"
#include <fstream>
#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
const unsigned DEFAULT_STRING_BUFFER = 100;
const unsigned DEFAULT_STEPS = 500000;
const double DEFAULT_STEP_SIZE = 0.001;
const unsigned DEFAULT_PRECISION = 5;
const unsigned DEFAULT_PARTICLES= 1;
const unsigned DEFAULT_MODE= 1;
const std::string DEFAULT_OUTPUT = "results";
const std::string DEFAULT_RESOURCES = "resources";
const unsigned DEFAULT_MAGPROF = 0;
const unsigned DEFAULT_NUM_POINTS = 10000;
const unsigned DEFAULT_PHI_ANGLE = 0;
const unsigned DEFAULT_PRINT_TYPE = 0;
const unsigned DEFAULT_DIMENSION = 1;


unsigned getPrintPrecisionFromArgs(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
        std::string tmp = argv[i];
        if(tmp == "-precision")
        {
            return static_cast<unsigned>(atoi(argv[i+1]));
        }
    }
    return DEFAULT_PRECISION;
}
unsigned getStepsFromArgs(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
        std::string tmp = argv[i];
        if(tmp == "-steps")
        {
            return static_cast<unsigned>(atoi(argv[i+1]));
        }
    }
    return DEFAULT_STEPS;
}
double getStepSizeFromArgs(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
        std::string tmp = argv[i];
        if(tmp == "-stepSize")
        {
            return strtod(argv[i+1], nullptr);
        }
    }
    return DEFAULT_STEP_SIZE;
}
void LoadParticles(const int& argc, char** argv, Coil& particles, const int length)
{
    bool found = false;
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
        std::string tmp = argv[i];
        if(tmp == "-particles")
        {
            loadFile(particles.x, particles.y, particles.z, length, argv[i+1]);
            found = true;
            break;
        }
    }
    if(!found)
    {
        printf("ERROR: particles path must be given!!\n");
        exit(1);
    }
}

std::string getResourcePath(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
        std::string param = argv[i];
        if("-resource" == param)
        {
            return std::string(argv[i+1]);
        }
    }
    return DEFAULT_RESOURCES;
}

unsigned getParticlesLengthFromArgs(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
        std::string tmp = argv[i];
        if(tmp == "-length")
        {
            return static_cast<unsigned>(atoi(argv[i+1]));
        }
    }
    printf("ERROR: length of particles path must be given!!\n");
    exit(1);
}
unsigned getModeFromArgs(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
        std::string tmp = argv[i];
        if(tmp == "-mode")
        {
            return static_cast<unsigned>(atoi(argv[i+1]));
        }
    }
    return DEFAULT_MODE;
}

std::string getJobId(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
        std::string tmp = argv[i];
        if(tmp == "-id")
        {
            return std::string(argv[i+1]);
        }
    }
    printf("ERROR: job id must be given!!\n");
    exit(1);
}

unsigned getMagneticProfileFromArgs(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
     	std::string tmp = argv[i];
        if(tmp == "-magnetic_prof")
        {
            return static_cast<unsigned>(atoi(argv[i+1]));
        }
    }
    return DEFAULT_MAGPROF;
}


unsigned getNumPointsFromArgs(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
     	std::string tmp = argv[i];
        if(tmp == "-magnetic_prof")
        {
            return static_cast<unsigned>(atoi(argv[i+2]));
        }
    }
    return DEFAULT_NUM_POINTS;
}

unsigned getAngleFromArgs(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
       	std::string tmp = argv[i];
        if(tmp == "-magnetic_prof")
        {
     	    return static_cast<unsigned>(atoi(argv[i+3]));
        }
    }
    return DEFAULT_PHI_ANGLE;
}

unsigned getDimension(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
     	std::string tmp = argv[i];
        if(tmp == "-magnetic_prof")
        {
            return static_cast<unsigned>(atoi(argv[i+4]));
        }
    }
    return DEFAULT_DIMENSION;
}

unsigned getPrintTypeFromArgs(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
     	std::string tmp = argv[i];
        if(tmp == "-print_type")
        {
            return static_cast<unsigned>(atoi(argv[i+1]));
        }
    }
    return DEFAULT_PRINT_TYPE;
}



int main(int argc, char** argv)
{


    std::string resourcePath; //Coil directory path
    unsigned steps; //Amount of simulation steps
    double stepSize; //Size of each simulation step
    
    /*Variables for magnetic profile diagnostic*/
    unsigned magprof; //Flag to control whether magnetic profile is computed
    unsigned num_points; //Number of sampling points for magnetic profile
    unsigned phi_angle; //Angle at which the magnetic profile will be computed
    /******************************************/
    
    unsigned precision; //TBD
    unsigned int length; //Amount of particles to simulate
    unsigned int mode; //Check divergence of simulation or not
    unsigned int print_type; //Printing flag: commas or tabs separators
    unsigned int dimension;
    std::string output; //Path of results directory
    std::string jobId; //JobID in the cluster
    std::ofstream handler;
    double startTime;
    double endTime;
	
    resourcePath = getResourcePath(argc, argv);
    steps = getStepsFromArgs(argc, argv);
    stepSize = getStepSizeFromArgs(argc, argv);
    precision = getPrintPrecisionFromArgs(argc, argv);
    length = getParticlesLengthFromArgs(argc, argv);
    mode = getModeFromArgs(argc, argv);
    print_type = getPrintTypeFromArgs(argc, argv);
    magprof = getMagneticProfileFromArgs(argc, argv);
    num_points = getNumPointsFromArgs(argc, argv);
    phi_angle = getAngleFromArgs(argc, argv);
    jobId = getJobId(argc, argv);

    dimension = getDimension(argc,argv);

    output = "results_" + jobId;
    createDirectoryIfNotExists(output);
    std::cout.precision(precision);
    std::cout << "Running with:" << std::endl;
    std::cout << "Resource Path=[" << resourcePath << "]." << std::endl;
    std::cout << "JobId=[" << jobId << "]." << std::endl;
    std::cout << "Steps=[" << steps << "]." << std::endl;
    std::cout << "Steps size=[" << stepSize << "]." << std::endl;
    std::cout << "Particles=[" << length << "]." << std::endl;
    std::cout << "Mode=[" << mode << "]." << std::endl;
    std::cout << "Output path=[" << output << "]." << std::endl;
    std::string file_name = "stdout_"+jobId+".log";
    
    handler.open(file_name.c_str());
    if(!handler.is_open()){
        std::cerr << "Unable to open stdout.log for appending. Nothing to do." << std::endl;
        exit(0);
    }

    handler << "Running with:" << std::endl;
    handler << "Steps=[" << steps << "]." << std::endl;
    handler << "Steps size=[" << stepSize << "]." << std::endl;
    handler << "Particles=[" << length << "]." << std::endl;
    handler << "Mode=[" << mode << "]." << std::endl;
    handler << "Output path=[" << output << "]." << std::endl;
	
    
    Coil particles;
    particles.x = static_cast<double*>(_mm_malloc(sizeof(double) * length, ALIGNMENT_SIZE));
    particles.y = static_cast<double*>(_mm_malloc(sizeof(double) * length, ALIGNMENT_SIZE));
    particles.z = static_cast<double*>(_mm_malloc(sizeof(double) * length, ALIGNMENT_SIZE));
    
    std::cout << "Loading Particles" << std::endl;
    LoadParticles(argc, argv, particles, length);
    std::cout << "Particles Loaded" << std::endl;


    const size_t sizeToAllocate = sizeof(double) * TOTAL_OF_GRADES_PADDED * TOTAL_OF_COILS;
    GlobalData data;
    data.coils.x = static_cast<double*>(_mm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    data.coils.y = static_cast<double*>(_mm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    data.coils.z = static_cast<double*>(_mm_malloc(sizeToAllocate, ALIGNMENT_SIZE));

    std::cout << "Loading Coils" << std::endl;
    load_coil_data(data.coils.x, data.coils.y, data.coils.z, resourcePath);
    std::cout << "Coil data loaded" << std::endl;
	
    data.e_roof.x = static_cast<double*>(_mm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    data.e_roof.y = static_cast<double*>(_mm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    data.e_roof.z = static_cast<double*>(_mm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    data.leng_segment = static_cast<double*>(_mm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    std::cout << "Preparing eroof" << std::endl;
    e_roof(data);
    std::cout << "e_roof data loaded" << std::endl;	    

    if(magprof != 0)
    {
	 std::cout << "Computing magnetic profiles" << std::endl;    
	 startTime = getCurrentTime();
         getMagneticProfile(data, num_points, phi_angle, output, dimension);
    	 endTime = getCurrentTime();
	 std::cout << "Total time computing magnetic profile=[" << (endTime - startTime) << "]." << std::endl;
    }

    startTime = getCurrentTime();
    std::cout << "Executing simulation" << std::endl;
    runParticles(data, output, particles, length, steps, stepSize, mode, print_type);
    endTime = getCurrentTime();
    std::cout << "Simulation finished" << std::endl;
    std::cout << "Total execution time=[" << (endTime - startTime) << "]." << std::endl;

    _mm_free(data.coils.x);
    _mm_free(data.coils.y);
    _mm_free(data.coils.z);
    _mm_free(data.e_roof.x);
    _mm_free(data.e_roof.y);
    _mm_free(data.e_roof.z);
    _mm_free(data.leng_segment);
    _mm_free(particles.x);
    _mm_free(particles.y);
    _mm_free(particles.z);

    handler << "Total execution time=[" << (endTime - startTime) << "]." << std::endl;
    handler.close();
    handler.open("stats.csv", std::ofstream::out | std::ofstream::app);
    if(!handler.is_open())
    {
        std::cerr << "Unable to open stats.csv for appending. Nothing to do." << std::endl;
        exit(0);
    }
    handler << jobId << "," << length << "," << steps << "," <<  stepSize << "," << output << "," << (endTime - startTime) << std::endl;
    handler.close();
    
    return 0;
}
