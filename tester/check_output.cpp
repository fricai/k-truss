#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <set>
using namespace std;

int main(int argc, char **argv){

    string org_file = argv[1];
    string out_file = argv[2];

    ifstream org_stream;
    ifstream out_stream;
    org_stream.open(org_file);
    out_stream.open(out_file);

    // comparing number of influencer

    int n_org;
    int n_out;



    string str_org;
    string str_out;

    
    getline(org_stream, str_org);
    n_org = stoi(str_org);
    map<int, set<int>> org_grps;
    for(int i = 0;i<2*n_org;i+=1){
	int curr;
	if(i%2==0){
	    getline(org_stream, str_org);
	    curr = stoi(str_org);
	    org_grps[curr] = set<int>();
	}
	else{
	    getline(org_stream, str_org);
	    istringstream iss(str_org);
	    
	    while(getline(iss, str_org, ' ')){
		int neb = stoi(str_org);
		org_grps[curr].insert(neb);
	    }
	    
	}
    }


    getline(out_stream, str_out);
    n_out = stoi(str_out);
    map<int, set<int>> out_grps;
    for(int i = 0;i<2*n_out;i+=1){
	int curr;
	if(i%2==0){
	    getline(out_stream, str_out);
	    curr = stoi(str_out);
	    out_grps[curr] = set<int>();
	}
	else{
	    getline(out_stream, str_out);
	    istringstream iss(str_out);
	    
	    while(getline(iss, str_out, ' ')){
		int neb = stoi(str_out);
		out_grps[curr].insert(neb);
	    }
	    
	}
    }
    bool matches = true;
    //

    if(n_out!=n_org){
	cout<<"DOES NOT MATCH"<<endl;
	return 0;
    }
    if(org_grps.size() != out_grps.size()){
	cout<<"DOES NOT MATCH"<<endl;
	return 0;
    }
    for(auto elem:org_grps){
	
	if(elem.second != out_grps[elem.first]){
	    cout<<"DOES NOT MATCH"<<endl;
	    return 0;
	}
	    
    }
    cout<<"MATCHES"<<endl;

    return 0;
}
