#include <iostream>
#include <thread>
#include <vector>
#include <algorithm>
#include <mutex>

#include <ql/quantlib.hpp>

using namespace QuantLib;

std::mutex mtx;
void generatorSobolSequenceSingle(SobolRsg& sobol,  std::vector<std::vector<double> >& results,std::vector<Time> times, Size stepRows, Size sampleCols)
{
    InverseCumulativeRsg<SobolRsg, InverseCumulativeNormal> gsg(sobol);
    BrownianBridge bridge(times);
    
    std::vector<double> increments(stepRows);
    for(Size j = 0; j< sampleCols ;++j)
    {
        const std::vector<double> sample = gsg.nextSequence().value;
        bridge.transform(sample.begin(), sample.end(), increments.begin());
        
        
        for(Size i = 1; i < stepRows;++i)
        {
            double stepValue = increments[i];
            mtx.lock();
            results[i][j] = stepValue;
            std::cout<< results[i][j] << std::endl;
            mtx.unlock();
        }
    }
};


int main(int argc, const char * argv[])
{
    //std::vector<Time> times = { 0.1, 0.2 , 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0};
    std::vector<Time> times_10y = { 1.0/3600.0, 2.0/3600.0, 3.0 / 3600.0, 4.0 /3600.0};
    
    std::vector<double> times = { 0.1, 0.2 , 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0};
    
    double constant = 1.0 /3600.0;
    
    std::transform(times.begin(), times.end(), times.begin(),
                       [ constant ](double element) { return element * constant; });
    
    TimeGrid grid(times.begin(), times.end());
    TimeGrid grid_10y(times_10y.begin(), times_10y.end());
    
    Size N = times.size();
    Size N_10y = times_10y.size();
    
    Size samples = 10;
    
    unsigned long  seed = 19910405;
    SobolRsg sobol(N, seed);
    SobolRsg sobol_10y(N_10y, seed);
    
    InverseCumulativeRsg<SobolRsg, InverseCumulativeNormal> gsg(sobol);
    InverseCumulativeRsg<SobolRsg, InverseCumulativeNormal> gsg_10y(sobol_10y);
    
    
    BrownianBridge bridge(times);
    BrownianBridge bridge_10y(times_10y);
    
    
    // temp : Z something //
    std::vector<double> temp(N);
    std::vector<double> temp_10y(N_10y);
    
    std::vector<double> temp_browinan(N);
    std::vector<double> temp_browinan_10y(N_10y);
    
    // paths //
    std::vector<std::vector<double> > vec_double(N, std::vector<double>(samples));
    std::vector<std::vector<double> > vec_dw(N, std::vector<double>(samples));
    
    Matrix mat_rn(N, samples);
    Matrix mat_dw(N, samples);
    
    SequenceStatistics stats1(N);
    
    
    for(Size j = 0; j < samples;++j)
    {
        const std::vector<double> sample = gsg.nextSequence().value;
        
        
        // temp : eta not Zi //
        bridge.transform(sample.begin(), sample.end(), temp.begin());
        
        // temp_brownian : Wi = Wi-1 + dWi, dWi = sqrt(dt) * eta  //
        bridge.transform(sample.begin(), sample.end(), temp_browinan.begin());
        
        temp_browinan[0] = temp[0] * std::sqrt(times[0]);
        std:: cout << j << "path(col) Browonian  " << std::endl;
        //time steps i//
        
        for(Size i =1 ; i <N ; ++i)
        {
            double dw =temp[i] * std::sqrt(times[i]-times[i-1]);
            temp_browinan[i] = temp_browinan[i-1] + dw;
            std::cout << i<< " step(row) dw" << dw << std::endl;
            std::cout << i<< " step(row) Wi " << temp_browinan[i] << std::endl;
            
            //std::cout << i << " step(row) dw " << dw << std::endl;
            vec_double[i][j] = temp_browinan[i];
            mat_rn(i,j) = temp_browinan[i];
            
            vec_dw[i][j] = dw;
            mat_dw(i,j) = dw;
            
        }
        //stats1.add(temp_browinan.begin(), temp_browinan.end());
    }
    
    std::cout << "vector || matrix for dw" << std::endl;
    for(int j = 0; j < samples; j++)
    {
        for(int i = 0; i < N;i++)
        {
            std::cout << j << "path , " << i << "step " << vec_dw[i][j] << " ||  (mat " << mat_dw(i,j) << ")"<< std::endl;
        }
    }
    
    
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    
    
    std::cout << "vector || matrix for single Wi" << std::endl;
    for(int j = 0; j < samples; j++)
    {
        for(int i = 0; i < N;i++)
        {
            std::cout << vec_double[i][j] << " || " << mat_rn(i,j) << std::endl;
        }
    }
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    
    std::cout << "testsuite ver " << std::endl;
    std::vector<double> testsuite(N);
    
    for(Size j = 0; j < samples;++j)
    {
        const std::vector<double> sample = gsg.nextSequence().value;
        
        bridge.transform(sample.begin(), sample.end(), testsuite.begin());
        
        
        std::cout << j<< " path(col) " << std::endl;
        
        testsuite[0] = testsuite[0] * std::sqrt(times[0]);
        //time steps//
        std:: cout << j <<  " Browonian " << std::endl;
        for(Size i =1 ; i <N ; ++i)
        {
            double dw =testsuite[i] * std::sqrt(times[i]-times[i-1]);
            testsuite[i] = testsuite[i-1] + dw;
            std::cout << i<< " step(row) " << testsuite[i] << std::endl;
            std::cout << i << " step(row) dw " << dw << std::endl;
        }
        //stats1.add(temp_browinan.begin(), temp_browinan.end());
    }
    
    
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    
    /*
    for(Size i = 0; i < samples;++i)
    {
        const std::vector<double> sample = gsg_10y.nextSequence().value;
        
        
        bridge_10y.transform(sample.begin(), sample.end(), temp_10y.begin());
        bridge_10y.transform(sample.begin(), sample.end(), temp_browinan_10y.begin());
        
        std::cout << i << " path(col) " << temp_10y[i] << std::endl;
        
        temp_browinan_10y[0] = temp_browinan_10y[0] * std::sqrt(times_10y[0]);
        //time steps//
        std:: cout << i <<  " Browonian " << std::endl;
        for(Size j =1 ; j <N_10y ; ++j)
        {
            temp_browinan_10y[j] = temp_browinan_10y[j-1] + temp_browinan_10y[j] * std::sqrt(times_10y[j]-times_10y[j-1]);
            std::cout << j << " step(row) " << temp_browinan_10y[j] << std::endl;
        }
        
    }
    */
    
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "multi threading " << std::endl;
    
    
    Size numThreads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    
    std::vector<std::vector<double> > temp_results(N, std::vector<double>(samples));
    
    std::cout << " test " << std::endl;
    std::vector<std::vector<double> > temp_results_test(N, std::vector<double>(samples));
    generatorSobolSequenceSingle(sobol, temp_results_test, times,N, samples);
    
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << " Mutlti Threading " << std::endl;
    int pathsPerThread = samples / numThreads;

    for (int i = 0; i < numThreads; ++i)
    {
        int start = i * pathsPerThread;
        int end = (i == numThreads - 1) ? samples : (i + 1) * pathsPerThread;
        threads.emplace_back(std::thread(generatorSobolSequenceSingle, std::ref(sobol), std::ref(temp_results),times, N, samples));
    }
    //for(int i = 0; i< numThreads; ++i)
    //{
    //    threads.push_back(std::thread(generatorSobolSequenceSingle, std::ref(sobol), std::ref(temp_results),times, N, samples));
    //}
    
    // Wait for all threads to finish
    std::for_each(threads.begin(), threads.end(), [](std::thread& t) { t.join(); });

    //for(auto& thread: threads)
    //{
    //    thread.join();
    //}
    
    for(Size i = 0; i < N; ++i)
    {
        for(Size j =0; j < samples; ++j)
        {
            std::cout << temp_results[i][j] << "\n";
        }
    }
    
    /*
    Size numThreads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    
    std::vector<std::vector<double> > temp_results(numThreads);
    for(Size i =0; i<numThreads;++i)
    {
        threads.push_back(std::thread(generatorSobolSequenceSingle, std::ref(sobol), std::ref(temp_results[i]), samples));
    }
    std::for_each(threads.begin(), threads.end(), [](std::thread& t) { t.join();});
    
    for(Size i = 0; i < numThreads; ++i)
    {
        std::cout << "Thread " << i << " : " << std::endl;
        for(Size row =0; row< N;++row)
        {
            
            std::cout << temp_results[][j] << " ";
        }
    }
     */
    return 0;
}
