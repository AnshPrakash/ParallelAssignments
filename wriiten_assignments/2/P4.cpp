#include <future>
#include <iostream>
#include <unistd.h>//for sleep function

using namespace std;


int Sum(future<int> &f){
	int sum=0;
	int n=f.get();
	for(int i=1;i<n+1;i++){
		sum+=i;
	}
	cout<<"Result of Thread "<<sum<<endl;
	cout<<"I have thread id "<<this_thread::get_id()<<endl;
	return(sum);
}

int main(int argc, char const *argv[])
{

	int x;
	promise<int> p;//Promises are container of future
	future<int> f=p.get_future();
 	future<int>	fu=async(launch::async,Sum,ref(f));//Sum function launch aynchronously in a separate thread(launch::async )
 	cout<<"Hello main thread here\n";
	cout<<"I have thread id "<<this_thread::get_id()<<"\n\n\n\n";
 	//Do some Work to give the value to Sum function
 	sleep(2);//doing some work 2s required (lot of work)
 	int val=100;//Ok Now I need the SUm of 1st 100 element
 	p.set_value(val);//As promised a value is given to p
 	//fu.get() will wait till the time the "Sum " exceution is complete,i.e,it is blocking
 	x=fu.get();//getting the returned value by thread2

	return 0;
}