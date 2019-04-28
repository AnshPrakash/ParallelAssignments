#include <bits/stdc++.h> 
#include "lab4_mpi.h"
#include <string.h> 
#include <malloc.h>
#include "mpi.h"
#include <iostream> 
#include <algorithm> 
#include <cmath>
#include <vector>
#include <set>

using namespace std;

void print_array(int* arr,int n){
	for(int i=0;i<n;i++){
		cout<<arr[i]<<" ";
	}
	cout<<"\n";
}


// Witness Function
int* Witness(char* Y,int p,int m){
	int pi = min(p,int(ceil(float(m)/2)));
	int* witness_array = (int*)calloc(pi,sizeof(int));
	witness_array[0] = 0;
	for(int i = 1;i<pi;i++){
		for(int k = 0;k<m-i;k++){
			if(Y[k]!=Y[i+k]){
				witness_array[i] = k;
				break;
			}
		}
	}
	return(witness_array);
}



int Duel(char* Z,int n,char* Y,int m,int* witness_array,int i,int j){
	int k = witness_array[j-i];
	if(j+k>=n ||  Z[j+k]!=Y[k]) return(i);
	else return(j);
}

bool pat_in_str(char* text,int n,char* pattern , int m,int pos){
	int j = 0;
	for(int i = pos;i<n,j<m ; i++){
		if(text[i]!=pattern[j]) return(false);
		j++;
	}
	return(true);
}

set<int> np_text_analysis(char* text,int n,char* pattern,int m,int* witness_array,int pi){
	//assumed pattern is non-periodic
	int blocks = int(ceil(n/ceil(float(m)/2)));
	int p = int(ceil(float(m)/2));
	int flag = 1;
	if (n%p == 0) flag =0;
	int i;
	int last;
	int potential_position[blocks];
	// cout<<"blocks "<<blocks<<"\n"; 
	for(int bi = 0;bi<blocks;bi++){
		i = (bi)*p;
		last = (bi + 1)*p;
		if(bi == blocks -2){
			last = (bi)*p + (1-flag)*p + flag*(n%p);
		}
		for(int j = i+1;j < last; j++){
			i = Duel(text,n,pattern,m,witness_array,i,j);
		}
		potential_position[bi] = i;
	}

	set<int> match_positions;
	for(int bi = 0 ; bi < blocks ; bi++){
		i = potential_position[bi];
		if(pat_in_str(text,n,pattern,m,i)){
			match_positions.insert(i);
			// cout<< i <<" ";
		}
	}
	// cout<<"\n";
	return(match_positions);
}

char* getsubstring(char* text, int n, int i, int j){
	//returns the substing starting at i and ending at j-1
	char* substr = (char*)(malloc(sizeof(char)*(j-i)));
	for(int k = 0 ; k < (j-i) ; k++){
		substr[k] = text[i+k];
	}
	return(substr);
}

char* generate(char* u, int l1, int k, char* v, int l2){
	char* gen = (char*)(malloc(sizeof(char)*(l1*k +l2)));
	for(int g =0 ;g<k;g++){
		for(int i =0;i<l1;i++){
			gen[(i+g*l1)] = u[i];
		}
	}
	for(int i =0;i<l2;i++){
		gen[(i+l1*k)] = v[i];
	}
	return(gen);
}


vector<int> check_cons_ones(vector<int> si,int num){
	vector<int> v(si.size());

	int count =0;
	// cout<<"num "<<num<<"\n";
	for (int i = si.size()-1; i >=0; --i){
		// cout<<"si" <<si[i]<<" ";
		if(si[i]==1) count++;
		else count = 0;
		// cout<<count<<"count"<<" ";
		if(count >= num ) v[i] = 1;
		else v[i] = 0;
	}
	// cout<<"\n";
	return(v);
}

bool somecondition(int j,int p){
	int i = j%p;
	if((i >=0 and i<p) ||(j<p)){
		return(true);
	}
	return(false);
}

vector<int> p_text_analysis(char* text,int n, char* pattern ,int m,int p){
	char* p_dash;
	p_dash = getsubstring(pattern,m,0,2*p-1);
	int* witness_array = Witness(p_dash,p,(2*p -1));

	int pi = p;//which is alway p = min(p,ceil(float(2*p -1)/2));
	set<int> pos = np_text_analysis(text,n,p_dash,2*p-1,witness_array,pi);
	char* u = getsubstring(pattern,n,0,p);
	
	

	int k = m/p;
	char* v = getsubstring(pattern,n,k*p,m);//till m-1

	


	char* u2v = generate(u,p,2,v,(m- k*p));
	int* M = (int*)(malloc(sizeof(int)*n));
	
	// for(int l =0;l<2*p+(m - k*p);l++){
		// cout<<u2v[l];
	// }
	// cout<<"\n";

	// cout<<"ps.count "<<pos.count(1)<<"\n";
	// cout<<"bool "<< pat_in_str(text,n,u2v,(2*p+(m-k*p)),1) <<"\n";
	for(int i = 0;i<n;i++){
		M[i] = 0;
		if(pos.count(i) && pat_in_str(text,n,u2v,(2*p+(m-k*p)),i) ){
			M[i] = 1;
		}
		// cout<<"M["<<i<<"]"<<M[i]<<" " ;
	}
	// cout<<"\n";
	// cout<<"k "<<k<<"\n";
	vector<vector<int>> C(p);
	for (int i = 0; i < p; ++i){
		vector<int> Si;
		int e = 0;
		while(i+e*p<n){
			Si.push_back(M[i+e*p]);
			e++;
		}
		C[i] = check_cons_ones(Si,k-1);
		
		// vector<int> ci(Si.size());
		// vector<int> conse_ones =check_cons_ones(Si,k-1) ;
		// for (int z  = 0; z < Si.size(); ++z){
		// 	cout<<"Consecutive ["<<z<<"]"<<conse_ones[z]<<" ";
		// }
		// cout<<"\n";
		// for (int z = 0; z <  Si.size() ; ++z){
		// 	cout<<"S["<<z<<"] "<<Si[z]<<" " ;
		// }
		// cout<<"\n";
		// for (int z = 0; z <  C[i].size() ; ++z){
		// 	cout<<"C["<<z<<"] "<<C[i][z]<<" " ;
		// }
		// cout<<"\n";
			
	}
	vector<int> MATCH;
	for(int j = 0;j <n-m +1;j++){
		if(somecondition(j,p)){
			int i = j%p;
			int l = j/p;
			if(C[i][l]) MATCH.push_back(j);
		}
	}
	return(MATCH);
}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void periodic_pattern_matching (int n, 
								char *text, 
								int num_patterns, 
								int *m_set, 
								int *p_set, 
								char **pattern_set,
								int **match_counts, 
								int **matches){
	
	int total_processes;
	MPI_Comm_size(MPI_COMM_WORLD, &total_processes);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	vector<vector<int>> local_matches;
	int start = (rank)*(num_patterns/(total_processes));
	int last = (rank+1)*(num_patterns/(total_processes));
	if (rank == total_processes-1) last = num_patterns;
	int local_size = 0;
	int local_pat_mat = 0;
	for(int i = start; i < last; i++){
		char* pattern;
		pattern =pattern_set[i];
		int p = p_set[i];
		local_size += 1;
		local_matches.push_back(p_text_analysis(text,n,pattern,m_set[i],p));
		local_pat_mat += local_matches[local_matches.size() -1].size();
		// cout<<"loc "<< local_matches[local_matches.size() -1].size()<<"\n";
	}

	
	int* data_in_each_process = (int*)(malloc(sizeof(int)*total_processes));
	int* matches_in_each_processes = (int*)(malloc(sizeof(int)*total_processes));
	(*match_counts)	= (int*)(malloc(sizeof(int)*num_patterns));

	int* local_match_count = (int*)(malloc(sizeof(int)*local_size));
	for (int i = 0; i < local_matches.size(); ++i){
		local_match_count[i] = local_matches[i].size();
	}
	
	// cout<<"my rank "<<rank<<" local_size "<<local_size<<"\n";
	MPI_Gather(&local_size, 1, MPI_INT, data_in_each_process , 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(&local_pat_mat, 1, MPI_INT,matches_in_each_processes, 1, MPI_INT, 0, MPI_COMM_WORLD);
	

	int* offsets;
	if(rank == 0){
		offsets = (int*)(calloc(total_processes,sizeof(int)));
		offsets[0] = 0;
		for (int i = 1; i < total_processes; ++i){
			offsets[i] = data_in_each_process[i-1] +offsets[i-1];
		}
	}

	MPI_Gatherv(local_match_count,local_size,MPI_INT
				,(*match_counts),data_in_each_process,offsets,MPI_INT,0,MPI_COMM_WORLD);

	int* local_data_matches = (int*)(malloc(sizeof(int)*local_pat_mat));
	int idx = 0;
	for (int i = 0; i < local_matches.size() ; ++i){
		for (int j= 0; j < local_matches[i].size(); ++j){
			// cout<<"local_pat_mat "<<local_pat_mat<<" idx "<<idx<<"\n";
			local_data_matches[idx] = local_matches[i][j];
			idx++;
		}
	}
	int* offsets2;
	int size_of_matches = 0;
	if(rank == 0){
		offsets2 = (int*)(calloc(total_processes,sizeof(int)));
		offsets2[0] = 0;
		for (int i = 1; i < total_processes; ++i){
			offsets2[i] = matches_in_each_processes[i-1] +offsets2[i-1];
		}
		for (int i = 0; i < total_processes; ++i){
			size_of_matches += matches_in_each_processes[i];
		}
			
	}

	// cout<<"local_pat_mat "<<local_pat_mat<<"\n";
	(*matches) =(int*)(malloc(sizeof(int)*size_of_matches));
	MPI_Gatherv(local_data_matches,local_pat_mat,MPI_INT
				,(*matches),matches_in_each_processes,offsets2,MPI_INT,0,MPI_COMM_WORLD);
	


	if (rank!=0){
		free(data_in_each_process);
		free(matches_in_each_processes);
		free((*match_counts));
		free((*matches));
	}

	// if (rank == 0){
		// *match_counts = (int*)(malloc(sizeof(int)*num_patterns));
		// (*matches) = (int*)(malloc(sizeof(int)*total));
		
		// vector<vector<int>> Match_bufer(num_patterns);
		
		//Checking the gathered result

		/////////////TESTING////////////////////////

		// for (int i = 0; i < total_processes; ++i){
		// 	cout<<data_in_each_process[i]<<" ";
		// }
		// cout<<"\n";
		// int count_1 =0;
		// for (int i = 0; i < num_patterns; ++i){
		// 	count_1 += (*match_counts)[i];
		// 	cout<<(*match_counts)[i]<<" ";
		// }
		// cout<<"\n";
		
		// for (int i = 0; i < size_of_matches; ++i){
		// 	cout<<(*matches)[i]<<" ";
		// }
		// cout<<"\n";
		// cout<<"Count 2 "<<size_of_matches;
		// cout<<"Count 1 "<< count_1;
		////////////TESTING OVER//////////////////




		// int* to_send_local_matches = (int*)(malloc(sizeof(int)*()));
		// for(int i = 0; i<local_matches.size(); i++){
		// 	m_count[i] = get<0>(local_matches[i]);
		// }

		//////////////////////////////////
		// int k = 0;
		// for (int i = 0; i < num_patterns; ++i){
		// 	for (int j = 0; j < Match_bufer[i].size(); ++j){
		// 		(*matches)[k] = Match_bufer[i][j];
		// 		k++;
		// 	}
		// }

		


	// }
	
	
	
}
