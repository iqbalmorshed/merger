#include<stdio.h>
#include<iostream>
#include<string>

using namespace std;


 void more_than_50p_remaining(string cigar);

int main(){

    //printf("%d %d",d,c);
    string str;
    //cout<<str.substr(2,50)<<end;
    while(1){
    cin>> str;
    more_than_50p_remaining(str);
    }


}

int total_matching(string cigar){

    int len_cigar = cigar.length();
    int sum_matching =0;
    int value =0;
    for(int i=0;i<len_cigar;i++){
        if(isdigit(cigar[i])){
                value = value*10 + (cigar[i] - '0');
                //printf("%d ",value);
        }
        else{
            if(cigar[i]=='M'){
                    sum_matching += value;

            }
             value =0;
        }
    }
    printf("sum of match = %d\N", sum_matching);

}


