#include <cmath>
#include <algorithm>
#include <cmath>

extern "C" {
	
	void zerov(double arr[],int a) {
		for(int i = 0; i < a; i++) {
			arr[i] = 0;
		}
		
	}
	
	double calculateMean(double arr[], int size) {
		double sum = 0.0;
		int count = 0; // Counter for non-zero elements
		for (int i = 0; i < size; ++i) {
			if (arr[i] != 0) { // Check if element is non-zero
				sum += arr[i];
				count++; // Increment count for non-zero elements
			}
		}
		if (count == 0) {
			return 0.0; // Return 0 if there are no non-zero elements
		} else {
			return sum / count; // Calculate mean using count of non-zero elements
		}
	}
	
	double calculateVariance(double arr[], int size) {
		double mean = calculateMean(arr, size);
		double variance = 0.0;
		for (int i = 0; i < size; ++i) {
			variance += (arr[i] - mean) * (arr[i] - mean);
		}
		return variance / size;
	}
	
	void calculateW(const double* S, double B, double* W, int size) {
		for (int i = 0; i < size; ++i) {
			W[i] = std::exp(-S[i] * B);
			}
		}
	
	double calculateFrost(const double* v, const double* W, int size) {
		double sum_v_times_W = 0.0;
		double sum_W = 0.0;
		for (int i = 0; i < size; ++i) {
			sum_v_times_W += v[i] * W[i];
			sum_W += W[i];
		}
		return sum_v_times_W / sum_W;
	}
	
	
	void frostfilter(double fdata[], const double data[], const double S[], const int prm[]) {
		int Z = prm[0];
		int X = prm[1];
		double az = prm[2]; 
		double ax = prm[3]; 
		
		int size = int(ax*az);
		int hax = ax/2;
		int haz = az/2;
		double z_pos;
		double x_pos;
		double v[size];
		double l_mean;
		double l_var;
		double B;
		double W[size];
		int count;
		for(int z = 0; z < Z; z++) {
			for(int x = 0; x < X; x++) {
				zerov(v,size);
				count = 0;
				for(int i = x-hax; i < x+hax; i++) {
					if(i>=0 && i<X){
						for(int j = z-haz; j < z+haz; j++) {
							if(j>=0 && j<Z){
								v[count] = data[j*X+i];
								count = count + 1;
							}
						}
					}
				}
				l_mean = calculateMean(v,size);
				l_var = calculateVariance(v,size);
				
				//B = Damp_fact * (l_var / (l_mean * l_mean));
				B = 1 * (l_var / (1+l_mean * l_mean));
				calculateW(S, B, W, size);
				fdata[z*X+x] = calculateFrost(v, W, size);
			}
		}
	}
}