///Bismillahir Rahmanir Rahim
#pragma warning(disable:4786)
#include <vector>
#include <list>
#include <map>
#include <set>
#include <deque>
#include <queue>
#include <stack>
#include <string>
#include <bitset>
#include <numeric>
#include <utility>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <functional>

#include <fstream>
#include <sstream>
#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cstring>
#include <cmath>
#include <cassert>
#include <ctime>
using namespace std;

#define MEMSET(dest,val) memset(dest,val,sizeof(dest))

#define STORE_COL 100
#define STORE_ROW 50
#define GRID_ROW 50

int match_limit = 60;
int merge_limit = 65;
void pre_processing();
void build_alignment_grid();
void update_grid_height(int row,int col);
void assemble();
void insert_align_grid(int first_pos,string cigar, string seq);
void print_align_grid();
void print_consensus(int startp, int endp);
void build_consensus();
void store_read(int point, char side, char read_store[STORE_ROW][STORE_COL]);
void build_store_profile(char read_store[STORE_ROW][STORE_COL], struct store_info store_profile[STORE_COL]);
void print_read_store(char read_store[STORE_ROW][STORE_COL]);
int find_right_merge_point(int left_point);
int find_distance_suggestion(int left_point,int right_point);
void merge_two_points(int left_point, int right_point, int distance_suggestion);
bool suggestion_merge(int suggestion, int left_point, int right_point);
bool shrink_merge(int suggestion, int left_point, int right_point);
bool expand_merge(int suggestion, int left_point, int right_point);
void print_assembly();
bool is_valid_for_expansion(int grid_col);
void print_high_insertion();
int total_matching(string cigar);


int read_length = 90;
int grid_height[1050] = {0};
int reference_length = 1000;
char left_read_store[STORE_ROW][STORE_COL] = {NULL};
char right_read_store[STORE_ROW][STORE_COL] = {NULL};
int new_grid_col = 0;

struct alignment_grid_info{ // alginment grid and colum starts from 1 (Not 0).
    char base; //could be 'A','T','G','C', 'D' (for deletion), NULL

    string str_insertion; //str_insertion is stored in the first position of the insertion. Ex- starting: 4 cigar: 2M2I3M will store the insertion at col 6
    string str_left_clip; //str_left_clip is stored in the first position of the aligned sequence. Ex- starting: 5 cigar: 3S2M will store the left clip at col 5
    string str_right_clip; //str_right_clip is stored in the last position of the aligned sequence. Ex- starting: 5 cigar: 3S2M4S will store the right clip at col 56



    bool is_insertion; // if this position/object holds str_insertion, then this variable will be true, otherwise false
    bool is_deletion;
    bool is_left_clip; // if this position/object holds str_left_clip, then this variable will be true, otherwise false
    bool is_right_clip; // if this position holds str_right_clip, then this variable will be true, otherwise false.
    bool is_clipping;
    bool is_starting_point;
    bool is_ending_point;

    alignment_grid_info(){
        base = NULL;
        str_insertion = "";
        is_insertion = false;
        is_deletion = false;
        is_left_clip = false;
        is_right_clip = false;
        is_clipping = false;
        is_starting_point = false;
        is_ending_point = false;
    }
}align_grid[GRID_ROW][1050]; //align grid row and col starts from 1. (not from 0)

struct base_count{
    char base;
    int cnt;
};

struct consensus_info{
    int cnt_A; // count of A
    int cnt_T;
    int cnt_G;
    int cnt_C;
    int cnt_D;
    int cnt_left_clip;
    int cnt_right_clip;
    int cnt_ins; // count of insertion
    int cnt_read;

    int perc_ins;
    int perc_del;
    int perc_left_clip;
    int perc_right_clip;

    char max_base;
    int cnt_max_base;
    int perc_max_base;

    consensus_info(){ //intialization
        cnt_A=0;cnt_T=0;cnt_G=0;cnt_C=0;cnt_left_clip=0;cnt_right_clip=0;cnt_left_clip=0;cnt_read=0;
        perc_ins =0;perc_del=0; perc_left_clip=0; perc_right_clip=0;
        max_base = NULL;
        cnt_max_base = 0;
        perc_max_base = 0;
    }

    void calculate_percentage_and_max_base(){

        if(cnt_read){
            perc_ins = (cnt_ins*100)/cnt_read;
            perc_left_clip = (cnt_left_clip*100)/cnt_read;
            perc_right_clip = (cnt_right_clip*100)/cnt_read;

            struct base_count bc_arr[5];
            bc_arr[0].base = 'A'; bc_arr[1].base = 'T'; bc_arr[2].base = 'G'; bc_arr[3].base = 'C'; bc_arr[4].base = 'D';
            bc_arr[0].cnt = cnt_A; bc_arr[1].cnt = cnt_T; bc_arr[2].cnt= cnt_G; bc_arr[3].cnt = cnt_C; bc_arr[4].cnt = cnt_D;


            for(int i=0;i<5;i++){
                if(bc_arr[i].cnt>cnt_max_base){
                    cnt_max_base = bc_arr[i].cnt;
                    max_base = bc_arr[i].base;
                }
            }

            perc_max_base = (cnt_max_base*100)/cnt_read;
        }


    }
}cons_arr[1050];

struct assembly_info{
    char base;
    int perc;
}assembled_unit;
vector <assembly_info> assembly;

struct store_info{
    int cnt_A;
    int cnt_T;
    int cnt_G;
    int cnt_C;
    int cnt_base;

    store_info(){
        cnt_A = cnt_T = cnt_G = cnt_C = cnt_base = 0;
    }
    void clear_store(){
        cnt_A = cnt_T = cnt_G = cnt_C = cnt_base = 0;
    }
}left_store_profile[STORE_COL], right_store_profile[STORE_COL];

struct merge_info{
    int cnt_A;
    int cnt_T;
    int cnt_G;
    int cnt_C;
    int cnt_base;
    int perc_max_base;
    char max_base;

    merge_info(){
        cnt_A = cnt_T = cnt_G = cnt_C = cnt_base =0;
        perc_max_base = 0;
        max_base = NULL;
    }

    void calculate_perc(){
        struct base_count bc_arr[4];
            bc_arr[0].base = 'A'; bc_arr[1].base = 'T'; bc_arr[2].base = 'G'; bc_arr[3].base = 'C';
            bc_arr[0].cnt = cnt_A; bc_arr[1].cnt = cnt_T; bc_arr[2].cnt= cnt_G; bc_arr[3].cnt = cnt_C;

            int cnt_max_base = 0;
            for(int i=0;i<4;i++){
                if(bc_arr[i].cnt>cnt_max_base){
                    cnt_max_base = bc_arr[i].cnt;
                    max_base = bc_arr[i].base;
                }
            }

            perc_max_base = (cnt_max_base*100)/cnt_base;

    }
};

int main(){

    //pre_processing();
    build_alignment_grid();
    //puts("hello");
    print_align_grid(); //for testing purpose
    build_consensus();
    //print_consensus(865,875);
    assemble();
    print_assembly();
    //store_read(647,'l',left_read_store);
    //print_read_store(left_read_store);
    //store_read(653,'r',right_read_store);
    //print_read_store(right_read_store);

    //printf("Distance suggestion: %d",find_distance_suggestion(908,921));
    //merge_two_points(908,921,23);
    //printf("new grid col: %d\n", new_grid_col);
    //print_read_store(left_read_store);
    //print_read_store(right_read_store);
    //print_high_insertion();

    return 0;
}
void print_high_insertion(){
    for(int j=1;j<reference_length;j++){
        if(cons_arr[j].perc_ins > 50){
            printf("insertion 50up in col %d, insertion perc %d\n", j, cons_arr[j].perc_ins);
        }
    }
}
void build_alignment_grid(){

    freopen("input.txt","r",stdin);
    freopen("output.txt","w",stdout);

    string str, header1, header2;

    getline(cin,header1);
    getline(cin,header2);

    int count =0,i,j;
    while(getline(cin,str)){
        stringstream ss;
        string qname, rname, cigar, rnext, seq,qual, quality_value_collection;
        int flag, pos, mapq, pnext, tlen;

        ss<<str;
        ss>>qname>>flag>>rname>>pos>>mapq>>cigar>>rnext>>pnext>>tlen>>seq>>qual;
        getline(ss,quality_value_collection);

        //cout<<qname<<" "<<pos<<" "<<cigar<<" "<<seq<<" "<<quality_value_collection<<endl;

        if(total_matching(cigar)< match_limit && pos>5){
            //cout<<"rejected cigar: "<<cigar<<endl;
        }
        else
            insert_align_grid(pos,cigar,seq);

    }



}
void insert_align_grid(int first_pos,string cigar, string seq){
    //cout<<first_pos<<" "<<cigar<<" "<<seq<<endl;
    int len_cigar = cigar.length();

    int value =0;
    int curr_grid_row = 1;
    int curr_grid_col = first_pos;
    int curr_seq_pos = 0;
    bool is_first_s =true;
    bool space_found = false;
    int i,j;

    for(i=1;i<50;i++){
        if(!align_grid[i][curr_grid_col].base){
            curr_grid_row=i;
            space_found = true;
            break;
        }
    }

    if(space_found){
        align_grid[i][curr_grid_col].is_starting_point= true;

        for(i=0;i<len_cigar;i++){
            if(isdigit(cigar[i])){
                value = value*10 + (cigar[i] - '0');
                //printf("%d ",value);
            }
            else{
                if(cigar[i]=='M'){

                    for( j=curr_grid_col;j<curr_grid_col+value; j++){
                        align_grid[curr_grid_row][j].base = seq[curr_seq_pos++];
                        //printf("%c ",align_grid[curr_grid_row][j].base);
                        update_grid_height(curr_grid_row,j);
                    }

                    curr_grid_col = j;
                    is_first_s=false;
                }
                else if(cigar[i]=='S'){
                    if(is_first_s){
                        align_grid[curr_grid_row][curr_grid_col].is_left_clip = true;
                        align_grid[curr_grid_row][curr_grid_col].is_clipping = true;
                        align_grid[curr_grid_row][curr_grid_col].str_left_clip = seq.substr(0,value);
                        is_first_s = false;
                        curr_seq_pos+=value;
                    }
                    else{
                        align_grid[curr_grid_row][curr_grid_col-1].is_right_clip = true;
                        align_grid[curr_grid_row][curr_grid_col-1].is_clipping = true;
                        //printf("hello1 %d %d",curr_seq_pos,read_length);
                        align_grid[curr_grid_row][curr_grid_col-1].str_right_clip = seq.substr(curr_seq_pos,value);
                        //printf("hello2");
                        curr_seq_pos+=value;
                    }
                }
                else if(cigar[i]=='D'){
                    for(j = curr_grid_col; j < curr_grid_col+value; j++){
                        align_grid[curr_grid_row][j].base = 'D';
                        align_grid[curr_grid_row][j].is_deletion = true;
                        update_grid_height(curr_grid_row,j);
                    }
                    curr_grid_col+=value;
                }
                else if(cigar[i]=='I'){
                    align_grid[curr_grid_row][curr_grid_col].is_insertion = true;
                    align_grid[curr_grid_row][curr_grid_col].str_insertion = seq.substr(curr_seq_pos,value);
                    curr_seq_pos+=value;
                    //doc: current_grid_col is not increased because next match will start from the current position.
                }
                value = 0;
            }

        }
        align_grid[curr_grid_row][curr_grid_col-1].is_ending_point = true;
    }
    //printf("curr seq: %d\n",curr_seq_pos);
    //if(curr_grid_row>align_grid_height[curr_grid_col])

}
void update_grid_height(int row,int col){
    if(row>grid_height[col])
        grid_height[col]=row;
}

void print_align_grid(){
    int i,j;
    //for(i=1;i<=270;i++)
        //printf("%d ",i);

    puts("");

    for(i=1;i<=49;i++){
        for(j=1;j<=500;j++){
            printf("%c ",align_grid[i][j].base);
        }
        puts("");
    }

    /*
    puts("");
    printf("printing grid height:\n");

    for(i=0;i<1000;i++){
        printf("grid col: %d height %d \n",i,grid_height[i]);
    }
    */
    //printf("special insertion: 28 650 %s\n",align_grid[28][650].str_insertion.c_str());
    printf("special left_clip: 16 121 %s\n",align_grid[16][121].str_left_clip.c_str());
}
void build_consensus(){

    int i,j;
    puts("building consensus....");
    for(j=1;j<reference_length;j++){
        for(i=1;i<=grid_height[j];i++){
            if(align_grid[i][j].base==NULL)
                continue;
            else if(align_grid[i][j].base=='A')
                cons_arr[j].cnt_A++;
            else if(align_grid[i][j].base=='T')
                cons_arr[j].cnt_T++;
            else if(align_grid[i][j].base=='G')
                cons_arr[j].cnt_G++;
            else if(align_grid[i][j].base=='C')
                cons_arr[j].cnt_C++;
            else if(align_grid[i][j].base=='D')
                cons_arr[j].cnt_D++;
            else{
                printf("Error in Align Grid. base is %c\n",align_grid[i][j].base);
            }

            if(align_grid[i][j].is_clipping){
                    //printf("inside_clipping col %d\n",j);
                if(align_grid[i][j].is_right_clip){
                    cons_arr[j+1].cnt_right_clip++; // Note that, count of right clip is stored in the count value of the next Consensus Array. Hence, j+1 is used. This is done in order to facilitate finding left merge point.
                    cons_arr[j+1].cnt_read++;
                }
                else
                    cons_arr[j].cnt_left_clip++;
            }

            if(align_grid[i][j].is_insertion)
                cons_arr[j].cnt_ins++;

            cons_arr[j].cnt_read++;
        }
        //printf("curr_col: %d\n",j);
        cons_arr[j].calculate_percentage_and_max_base();

    }



}
void print_consensus(int startp, int endp){
    for(int j=startp;j<=endp;j++){
        printf("Col: %d\n",j);
        printf("max base: %c %d\n",cons_arr[j].max_base, cons_arr[j].perc_max_base);
        printf("right clip perc & conut: %d %d\n",cons_arr[j].perc_right_clip, cons_arr[j].cnt_right_clip);
        printf("left clip perc: %d\n",cons_arr[j].perc_left_clip);
        printf("insertion perc: %d\n",cons_arr[j].perc_ins);
        printf("# of reads: %d, # of insertion: %d, # of left clip %d, # of right clip %d\n",cons_arr[j].cnt_read, cons_arr[j].cnt_ins, cons_arr[j].cnt_left_clip, cons_arr[j].cnt_right_clip);
        puts("");
    }
}

void assemble(){

    int grid_col = 1;
    //bool is_found_left_merge_point = false;

    while(grid_col<= reference_length){

        printf("current grid col %d\n",grid_col);

        if(cons_arr[grid_col].perc_ins >=75){
            puts("insertion possibility....................................................");
            for(int i=1;i<grid_height[grid_col];i++){
                if(align_grid[i][grid_col].base && align_grid[i][grid_col].is_insertion){
                    printf("Insertion occured in col %d ============== \n", grid_col);
                    int len = align_grid[i][grid_col].str_insertion.length();

                    for(i=0;i<len;i++){
                        assembled_unit.base = align_grid[i][grid_col].str_insertion[i];
                        assembled_unit.perc = 0;
                        assembly.push_back(assembled_unit);
                    }
                    break;
                }
            }
        }

        if(cons_arr[grid_col].perc_max_base >= 75 ){
            assembled_unit.base = cons_arr[grid_col].max_base;
            assembled_unit.perc = cons_arr[grid_col].perc_max_base;

            assembly.push_back(assembled_unit);
            grid_col++;

        }
        else if(is_valid_for_expansion(grid_col) ){
            // This are needs work. Program won't enter here because max base here will always be greater.
            printf("Entered for expansion======");
            //merge_two_points(grid_col-1,grid_col, 0);
                //printf("Expand merge occured due to left clipping at col %d============",grid_col);
            //grid_col = new_grid_col;

            grid_col++; //temp
        }
        else{
            int left_merge_point, right_merge_point, distance_suggestion;

            left_merge_point = grid_col -1;
            right_merge_point = find_right_merge_point(left_merge_point);

            if(right_merge_point >= reference_length){
                    puts("right merge point out of range");
                    break;
            }
            distance_suggestion = find_distance_suggestion(left_merge_point,right_merge_point);
            merge_two_points(left_merge_point, right_merge_point, distance_suggestion);

            if(new_grid_col< grid_col){
                printf("lesser: new grid %d, curr grid: %d left_point %d right point %d suggestion: %d\n",new_grid_col,grid_col, left_merge_point, right_merge_point, distance_suggestion);
                print_read_store(left_read_store);
                print_read_store(right_read_store);
                break;
            }
            grid_col = new_grid_col;
        }
    }
}
bool is_valid_for_expansion(int grid_col){
    if(grid_col > 5 && (cons_arr[grid_col].perc_left_clip > 15) )
        return true;
}
int find_right_merge_point(int left_point){

    int col;
    for(col = left_point+1; col<=reference_length ;col++){

        if(cons_arr[col].perc_max_base>=85 && cons_arr[col].cnt_max_base>= cons_arr[left_point].cnt_max_base )
            break;
    }

    return col;

}
int find_distance_suggestion(int left_point,int right_point){

    int i,j;
    vector <int> suggestion(GRID_ROW, -1); // -1 means No Idea

    for(i= 1;i < grid_height[left_point];i++){
        if(align_grid[i][left_point].base && !align_grid[i][left_point].is_right_clip){
            suggestion[i]=0;
            for( j=left_point+1;j<right_point;j++){

                if(align_grid[i][j].is_right_clip){
                    suggestion[i] = -1; //this row has no idea.
                    break;
                }
                if(align_grid[i][j].base != 'D')
                        suggestion[i]++;

                if(align_grid[i][j].is_insertion)
                        suggestion[i] += align_grid[i][j].str_insertion.length();

            }
        }
    }

    //for(i=1;i< grid_height[left_point];i++)
        //printf("%d ",suggestion[i]);

    map <int,int> freq_mp;
    map <int,int>::iterator mit;

    for(i=1; i < grid_height[left_point]; i++){
        if(suggestion[i]!= -1){
            freq_mp[suggestion[i]]++;
        }
    }
    int max_freq =0;
    int max_suggestion = 0;
    for(mit= freq_mp.begin();mit!=freq_mp.end();mit++){
        if(mit->second > max_freq){
            max_freq = mit->second;
            max_suggestion = mit->first;
        }
    }
    return max_suggestion;
}
void merge_two_points(int left_point, int right_point, int distance_suggestion){

    store_read(left_point, 'l', left_read_store);
    store_read(right_point, 'r', right_read_store);
    build_store_profile(left_read_store, left_store_profile);
    build_store_profile(right_read_store, right_store_profile);

    printf("Merge happening =====> left point: %d, right point: %d suggestion: %d\n",left_point, right_point, distance_suggestion);

    if(suggestion_merge(distance_suggestion, left_point, right_point)){
        printf("Suggestion merge worked for sug: %d left: %d right: %d\n",distance_suggestion, left_point, right_point);
    }
    else if(shrink_merge(distance_suggestion, left_point, right_point)){
        printf("Shrink merge worked for sug: %d left: %d right: %d\n",distance_suggestion, left_point, right_point);
    }
    else if(expand_merge(distance_suggestion, left_point, right_point)){
        printf("Expand merge worked for sug: %d left: %d right: %d\n",distance_suggestion, left_point, right_point);
    }
    else{
        printf("All merger failed\n");
        new_grid_col = right_point;
    }

}

void store_read(int point, char side, char read_store[STORE_ROW][STORE_COL]){

    int left_starting, right_starting;
    memset(read_store, 0, sizeof(read_store[0][0]) * STORE_ROW * STORE_COL);

    if(side == 'l'){
        left_starting = point-1;
        right_starting = point;
    }
    else{
        left_starting = point;
        right_starting = point+1;
    }

    int char_limit = STORE_COL/2; //char limit on each side;
    int store_row = 0;

    for(int grid_row = 1; grid_row <= grid_height[point]; grid_row++){
        int grid_col = left_starting;
        if(align_grid[grid_row][point].base){

            //for left part of merge point
            int cnt_char =0;
            int store_col = char_limit-1;
            //if(side == 'r')printf("store row: %d store col: %d grid row: %d grid col %d point %d\n",store_row,store_col,grid_row,grid_col,point);
            while(cnt_char < char_limit){
                if(align_grid[grid_row][grid_col].base != 'D'){
                    read_store[store_row][store_col]= align_grid[grid_row][grid_col].base;
                    cnt_char++;
                    store_col--;

                    if(align_grid[grid_row][grid_col].is_insertion){
                        int len = align_grid[grid_row][grid_col].str_insertion.length();
                        while(cnt_char < char_limit && len>0){
                            read_store[store_row][store_col]= align_grid[grid_row][grid_col].str_insertion[len-1];
                            cnt_char++;
                            store_col--;
                            len--;
                        }
                    }

                    if(align_grid[grid_row][grid_col].is_left_clip){
                        int len = align_grid[grid_row][grid_col].str_left_clip.length();
                        while(cnt_char< char_limit && len>0){
                            read_store[store_row][store_col]= align_grid[grid_row][grid_col].str_left_clip[len-1];
                            cnt_char++;
                            store_col--;
                            len--;
                        }
                        break;
                    }

                    if(align_grid[grid_row][grid_col].is_starting_point)
                        break;
                }
                grid_col--;
            }

            //for right part of merge point
            cnt_char =0;
            grid_col = right_starting;
            store_col = char_limit;
            while(cnt_char < char_limit){
                if(align_grid[grid_row][grid_col].base!='D'){

                    if(align_grid[grid_row][grid_col].is_insertion){
                        int len = align_grid[grid_row][grid_col].str_insertion.length();
                        int i=0;
                        while(cnt_char< char_limit && i<len){
                            read_store[store_row][store_col]= align_grid[grid_row][grid_col].str_insertion[i];
                            cnt_char++;
                            store_col++;
                            i++;
                        }

                    }

                    if(cnt_char < char_limit){
                        read_store[store_row][store_col]= align_grid[grid_row][grid_col].base;
                        cnt_char++;
                        store_col++;

                    }

                    if(align_grid[grid_row][grid_col].is_right_clip){       //is_right_clip used becasue left clip is not possible in this part.
                        int len = align_grid[grid_row][grid_col].str_right_clip.length();
                        int i=0;
                        while(cnt_char< char_limit && i<len){
                            read_store[store_row][store_col]= align_grid[grid_row][grid_col].str_right_clip[i];
                            cnt_char++;
                            store_col++;
                            i++;
                        }

                        break;
                    }

                    if(align_grid[grid_row][grid_col].is_ending_point)
                        break;
                }

                grid_col++;
            }

            ///////// completed one store row
            store_row++;

        }
    }
}
void print_read_store(char read_store[STORE_ROW][STORE_COL]){
    int i,j;
    for(i=0;i<30;i++){
        for(j=0;j<STORE_COL;j++){
            printf("%c",read_store[i][j]);
            //if(left_read_store[i][j]==NULL)
                //printf("%c",'x');
        }
        puts("");
    }
}

void build_store_profile(char read_store[STORE_ROW][STORE_COL], struct store_info store_profile[STORE_COL]){

    for(int j=0;j<STORE_COL;j++)
        store_profile[j].clear_store();

    int i,j;
    for(j=0;j<STORE_COL;j++){
        for(i=0;i<STORE_ROW;i++){
            if(read_store[i][j]){
                if(read_store[i][j] == 'A')
                    store_profile[j].cnt_A++;
                else if(read_store[i][j] == 'T')
                    store_profile[j].cnt_T++;
                else if(read_store[i][j] == 'G')
                    store_profile[j].cnt_G++;
                else
                    store_profile[j].cnt_C++;

                store_profile[j].cnt_base++;
            }
        }
    }
}
bool suggestion_merge(int suggestion, int left_point, int right_point){

    int allowence = 1;
    int start_point, end_point, right_store_col;
    if(suggestion <=0){ //0 means no distance. negative means there's overlap. positive means there is distance.
        start_point = STORE_COL/2 - 4 + suggestion;  //start from 4 bp left from left point to 4 bp right to right  point. So if
        end_point = STORE_COL/2 +5 + suggestion;
    }
    else{
        start_point = STORE_COL/2 -4;
        end_point = STORE_COL/2 +5 + suggestion;
    }

    if(start_point < 0){
        start_point =0;
        printf("Left limit exceed");
    }

    if(end_point >= STORE_COL){
        end_point = STORE_COL -1;
        printf("Right limit exceed");
    }

    struct merge_info merge_arr[STORE_COL];

    //printf("start_point: %d end_point: %d\n",start_point, end_point);
    int merge_col = 0, j;
    for(j=start_point;j<=end_point;j++){
        right_store_col = j - 2 - suggestion;

        merge_arr[merge_col].cnt_A = left_store_profile[j].cnt_A + right_store_profile[right_store_col].cnt_A;
        merge_arr[merge_col].cnt_T = left_store_profile[j].cnt_T + right_store_profile[right_store_col].cnt_T;
        merge_arr[merge_col].cnt_G = left_store_profile[j].cnt_G + right_store_profile[right_store_col].cnt_G;
        merge_arr[merge_col].cnt_C = left_store_profile[j].cnt_C + right_store_profile[right_store_col].cnt_C;
        merge_arr[merge_col].cnt_base = left_store_profile[j].cnt_base + right_store_profile[right_store_col].cnt_base;


        merge_arr[merge_col].calculate_perc();
        printf("merge col %d rotating\n left_store col: %d, right store col: %d suggestion: %d perc max base %d max_base %c\n",merge_col,j,j- 2 - suggestion, suggestion,merge_arr[merge_col].perc_max_base,merge_arr[merge_col].max_base);

        if(merge_arr[merge_col].perc_max_base < merge_limit){
            if(!allowence)
                return false;
            else
                allowence--;
        }

        merge_col++;

    }

    if(suggestion == 0){
        new_grid_col = right_point;
    }
    else if(suggestion < 0){
        new_grid_col = right_point - suggestion;
    }
    else{
        new_grid_col = right_point;
        for(int i= 5; i < 5+suggestion; i++){
            assembled_unit.base = merge_arr[i].max_base;
            assembled_unit.perc = merge_arr[i].perc_max_base;
            assembly.push_back(assembled_unit);
        }
    }

    return true;
}
bool shrink_merge(int suggestion, int left_point, int right_point){

    for(int i = suggestion-1; i> -10;i--){
        if(suggestion_merge(i, left_point, right_point))
            return true;
    }
    return false;
}
bool expand_merge(int suggestion, int left_point, int right_point){

    for(int i = suggestion+1; i < 10;i++){
        if(suggestion_merge(i, left_point, right_point))
            return true;
    }
    return false;
}
void print_assembly(){


    int len = assembly.size();
    printf("Assembled String: Length: %d\n", len);
    for(int i=0;i<len;i++){
        if(assembly[i].base != 'D')
            printf("%c",assembly[i].base);
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
    return sum_matching;
    //printf("sum of match = %d\N", sum_matching);

}
