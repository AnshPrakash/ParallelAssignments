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

	

	*match_counts = (int*)(malloc(sizeof(int)*num_patterns));
	int total = 0;
	vector<vector<int>> Match_bufer(num_patterns);
	for(int i = 0;i<num_patterns;i++){
		char* pattern =(char*)(malloc(sizeof(char)*m_set[i])) ;
		int p = p_set[i];
		for(int j =0;j<m_set[i];j++){
			pattern[j] = (*pattern_set)[i*m_set[i]+j];
		}
		Match_bufer[i] =  p_text_analysis(text,n,pattern,m_set[i],p);
		(*match_counts)[i] = Match_bufer[i].size();
		total += Match_bufer[i].size();
		free(pattern);
	}
	(*matches) = (int*)(malloc(sizeof(int)*total));
	int k = 0;
	for (int i = 0; i < num_patterns; ++i){
		for (int j = 0; j < Match_bufer[i].size(); ++j){
			(*matches)[k] = Match_bufer[i][j];
			k++;
		}
	}
}
