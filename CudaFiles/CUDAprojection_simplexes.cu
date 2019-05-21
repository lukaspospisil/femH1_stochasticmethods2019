__device__ void CUDAprojection_simplexes_sort_quickSort_partition(int *pi, double *x, int low, int high, int t, int T)
{
    /* pivot (Element to be placed at right position) */
    double pivot = x[high*T + t];
 
    int i = (low - 1);  /* Index of smaller element */
	double swap;

    for(int j = low; j <= high - 1; j++)
    {
        // If current element is smaller than or
        // equal to pivot
        if(x[j*T + t] <= pivot)
        {
            i++;    // increment index of smaller element
            
			swap = x[i*T+t];
			x[i*T+t] = x[j*T+t];
			x[j*T+t] = swap;
        }
    }

    swap = x[(i+1)*T+t];
	x[(i+1)*T+t] = x[high*T+t];
	x[high*T+t] = swap;

    *pi = i + 1;
}

__device__ void CUDAprojection_simplexes_sort_quickSort(double *x, int low, int high, int t, int T)
{
    if (low < high)
    {
        /* pi is partitioning index, arr[pi] is now
           at right place */
		int pi;
		CUDAprojection_simplexes_sort_quickSort_partition(&pi, x, low, high, t,T);

		CUDAprojection_simplexes_sort_quickSort(x, low, pi - 1, t,T);  /* Before pi */
		CUDAprojection_simplexes_sort_quickSort(x, pi + 1, high, t,T); /* After pi */
    }
}

__device__ void device_sort_bubble(double *x_sorted, int t, int T, int K){
	int i;
	int m=K;
	int mnew;
	double swap;

	while(m > 0){
		/* Iterate through x */
		mnew = 0;
		for(i=1;i<m;i++){
			/* Swap elements in wrong order */
			if (x_sorted[i*T+t] < x_sorted[(i - 1)*T + t]){
				swap = x_sorted[i*T + t];
				x_sorted[i*T + t] = x_sorted[(i - 1)*T + t];
				x_sorted[(i - 1)*T + t] = swap;
				mnew = i;
			}
        }
		m = mnew;
	}
}

__global__ void CUDAprojection_simplexes( double *X, 
                   double *Y,
                   int T, int K ) {

	int t = blockIdx.x*blockDim.x + threadIdx.x;
	
	int k;

	if(t<T){
		bool is_inside = true;
		double sum = 0.0;
	
		/* control inequality constraints */
		for(k = 0; k < K; k++){ // TODO: could be performed parallely  
			if(X[k*T+t] < 0.0){
				is_inside = false;
			}
			sum += X[k*T + t];
			
			Y[k*T + t] = X[k*T + t];
		}

		/* control equality constraints */
		if(sum != 1){ 
			is_inside = false;
		}

		/* if given point is not inside the feasible domain, then do projection */
		if(!is_inside){
			int j,i;
			/* compute sorted x_sub */
			double sum_y;

//			CUDAprojection_simplexes_sort_bubble(Y,t,T,K);
            CUDAprojection_simplexes_sort_quickSort(Y, 0, K-1, t, T);

			/* now perform analytical solution of projection problem */	
			double t_hat = 0.0;
			i = K - 1;
			double ti;

			while(i >= 1){
				/* compute sum(y) */
				sum_y = 0.0;
				for(j=i;j<K;j++){ /* sum(y(i,n-1)) */
					sum_y += Y[j*T + t];
				}
				
				ti = (sum_y - 1.0)/(double)(K-i);
				if(ti >= Y[(i-1)*T + t]){
					t_hat = ti;
					i = -1; /* break */
				} else {
					i = i - 1;
				}
			}

			if(i == 0){
				t_hat = (sum-1.0)/(double)K; /* uses sum=sum(x_sub) */
			}
    
			for(k = 0; k < K; k++){ // TODO: could be performed parallely  
				/* (*x_sub)(i) = max(*x_sub-t_hat,0); */
				ti = X[k*T + t] - t_hat;	
				if(ti > 0.0){
					X[k*T + t] = ti;
				} else {
					X[k*T + t] = 0.0;
				}
			}
		}
		
	}

	/* if t >= T then relax and do nothing */
}
