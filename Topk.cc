//
//  TopK.cc
//
//
//  Created by Zichuan Zou on 3/16/14.
//
// 

#include <vector>
#include <cctype>
#include <string>
#include <iomanip>
#include <stdint.h>
#include <fstream>
#include <iterator>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <sys/time.h>
#include <unordered_map>

using namespace std;
//using namespace tr1;

#define UNDEFROW -100000

struct Node;

//initialization of globle variable
unsigned maxLen = 0;
unsigned minLen = 0x7FFFFFFF;

vector<string> data;  //the vector for the collection of strings
vector<int> prefixLen;  //the prefixs of each string
vector< pair<unsigned, unsigned> > res;    //result set, .fisrt is the no. of string after sort; the .second is the edit distance of the string
vector< pair<string, unsigned> > queryInfos;     //the query and its k value
unordered_map<char, vector<Node> > * charMapId;    //in the trie, I have each char as a node of the trie with its corresponding information contains in the Node struture


//node on trie. this is a little bit different from paper.
struct Node {
  pair<unsigned, unsigned> idRange;
  pair<unsigned, unsigned> lengths;
  unordered_map<char, int> children;
  Node(unsigned lo_id, unsigned up_id, unsigned lo_len, unsigned up_len) {
    idRange = make_pair(lo_id, up_id);
    lengths = make_pair(lo_len, up_len);
  }
};

//endpoint: the location of last pivatol entries. contains the range of ID of endpoint. notice there are 2 different ways to initialize it
struct EndPoint {
  int row;
  pair<unsigned, unsigned> range;
  EndPoint(int _row, unsigned lo, unsigned up)
    : row(_row) { range = make_pair(lo, up); }
  EndPoint() {
    row = UNDEFROW;
    range = make_pair(0, data.size() - 1);
  }
};

//declaration
void TopkSearch(const string & query, const unsigned topk);

//compare 2 strings, return true if s2>s1
bool comp(const string & s1, const string & s2) { return s1 < s2; }

//below is main function
/////////////////////////////////////////////////////////////////////////////


int main(int argc, char ** argv) {
  char * dataFilename = argv[1];
  char * queryFilename = argv[2];

    
    
  string line;
    
  //read the datafile
  ifstream dataFile(dataFilename, ios::in);
  while (getline(dataFile, line, '\n') != NULL) {
    if (maxLen < line.length()) maxLen = line.length();//thus I know the Max and Min length of these strings
    if (minLen > line.length()) minLen = line.length();
    data.push_back(line); //push them into data vector
  }
  dataFile.close();

  //read the queryfile, this can input multiple queries at the same time
  ifstream queryFile(queryFilename, ios::in);
  while (getline(queryFile, line, '\n') != NULL) {
    string::size_type pos = line.rfind(' ');
    if (pos == string::npos) continue;
    string query = line.substr(0, pos);
    unsigned topK = atoi(line.substr(pos + 1).c_str());//return the k of topk
    queryInfos.push_back( pair<string,unsigned>(query, topK) );
  }
  queryFile.close();

    
  //sort all the input records so that the shortest string will be at fisrt place
  sort(data.begin(), data.end(), comp);

    
    
  // create the index, its size is the max length of strings. this is the trie
  charMapId = new unordered_map<char, vector<Node> > [maxLen];

  //vector<int> prefixLen; push the 0 to prefixLen. this initialization means there's no differences between 2 strings in string set so far
  prefixLen.push_back(0);
    
    //for the 1st string in set of Strings
  for (int lp = 0; lp < data[0].length(); lp++) {
    charMapId[lp][data[0][lp]].push_back(Node(0, 0, data[0].length(), data[0].length()));//Node(unsigned lo_id, unsigned up_id(idRange), unsigned lo_len, unsigned up_len(length))
      if (lp > 0) charMapId[lp-1][data[0][lp-1]][0].children.insert(unordered_map<char, int>::value_type(data[0][lp], 0));
  } // unordered_map<char, int> children; put the data[0][lp] into into this node's children。chars after the 1st char will be its children
                                                                                                       
    //traverse each string
  for (int iter = 1; iter < data.size(); iter++) {
    bool flag = true;//if the string is different from the last string, then it should go the other way.
    int len = data[iter].length();//len = length of each string
      
    for (int lp = 0; lp < len; lp++) {
      char ch = data[iter][lp];// iter: the string that this char resides on; lp: it's position index
      if (flag == true) {
        if (data[iter][lp] != data[iter - 1][lp]) {//if this string is different in chars in the same place with its previous string
          prefixLen.push_back(lp);//where is this difference located in strings
          flag = false; //false mean there's a difference exists
          charMapId[lp][ch].push_back(Node(iter, iter, len, len));//add this different char in the same location with previous char
            //add the children for previous char
          if (lp > 0) {
            int lastOne = charMapId[lp][ch].size() - 1;
            charMapId[lp-1][data[iter][lp-1]].back().children.insert(make_pair(ch, lastOne));
          }
        }
        //if this string is same the former string in char, then do
        else {
          charMapId[lp][ch].back().idRange.second = iter;
          if (charMapId[lp][ch].back().lengths.first > len) 
            charMapId[lp][ch].back().lengths.first = len;
          if (charMapId[lp][ch].back().lengths.second < len) 
            charMapId[lp][ch].back().lengths.second = len;
        }
      }
      
      //if the flag is false
      else {
        charMapId[lp][ch].push_back(Node(iter, iter, len, len));
        int lastOne = charMapId[lp][ch].size() - 1;
        charMapId[lp-1][data[iter][lp-1]].back().children.insert(make_pair(ch, lastOne));
      }
    }
  }
                                                                                                       
  //set the output file if user didn't set
  ofstream output(argc > 3 ? argv[3] : "output.txt");
  
  //output the time and query results
  timeval tbegin, tend;
  gettimeofday(&tbegin, NULL);
                                                                                                       
  //vector< pair<string, unsigned> > queryInfos; vector< pair<unsigned, unsigned> > res
  for (int i = 0; i < queryInfos.size(); i++) {
    //each iteration should clean the res so that the result of next query could be stored
      res.clear();
    printf("======== QUERY: %s %d ======\n", queryInfos[i].first.c_str(), queryInfos[i].second);
      //do the search for every query
    TopkSearch(queryInfos[i].first, queryInfos[i].second);
      //read the string in the res one by one and print it
    for (int j = 0; j < res.size(); j++)
      printf("%3d : %d %s %d\n", j + 1, res[j].first, data[res[j].first].c_str(), res[j].second);
    printf("===================== END OF QUERY ===================\n");
    //print the ed of last element in res
      output << res.back().second << endl;
  }
                                                                                                       
                                                                                                       
  gettimeofday(&tend, NULL);
    
  //output the elapsd time
  printf("# Time Elapsd : %f\n", tend.tv_sec - tbegin.tv_sec + (tend.tv_usec - tbegin.tv_usec) * 1.0 / CLOCKS_PER_SEC);
  output << tend.tv_sec - tbegin.tv_sec + (tend.tv_usec - tbegin.tv_usec) * 1.0 / CLOCKS_PER_SEC << endl;
  
    //clean the buffuer and temp data
  output.close();
  delete [] charMapId;
  return 0;
}


////////////////////////////////////////////////////////////////////////////
//main function over


//return the max one among 3 numbers
inline int maximum(int a, int b, int c) { return a >= b ? (a >= c ? a : c) : (b >= c ? b : c); }

//if node p1's largest range number is smaller than p2
bool CompareLower(const Node & p1, const unsigned p2) { return p1.idRange.second < p2; }

//if p2's smallest number is larger than p1
bool CompareUpper(const unsigned p1, const Node & p2) { return p1 < p2.idRange.first; }

//RollDown(diagonals[!f][id + 1], dia, row, lo, up, query, ed, topk): narrow down the scope of range and row
bool RollDown(vector<EndPoint> & toFillIn, int k, int row, unsigned low, unsigned upp, const string & query, int ed, int topk) {
  unsigned lo = low;
  unsigned up = upp;
  vector<int> record;
  unordered_map<char, int>::iterator mit;
  vector<Node>::iterator begin, it;

  // first item.idRange.second >= lo
  begin = lower_bound(charMapId[row+1][query[row+1+k]].begin(), charMapId[row+1][query[row+1+k]].end(), lo, CompareLower);

  // CHECK THE CORRECTNESS. while the low bound is not the same with lo, we need to deduct it
  while (begin != charMapId[row+1][query[row+1+k]].end() && begin->idRange.first <= upp) {
    if (begin->idRange.first > lo) {
      toFillIn.push_back(EndPoint(row, lo, begin->idRange.first - 1));
      lo = begin->idRange.first;
    }//choose the lowest upper bound
    if (begin->idRange.second <= upp) {
      record.push_back(begin->idRange.second);
      up = begin->idRange.second;
    } else {
      record.push_back(upp);
      up = upp;
    }
    ++row;
    it = begin;
    
    while (row + 1 + k < query.length()) {
      mit = it->children.find(query[row+1+k]);
      if (mit == it->children.end()) break;
      it = charMapId[row+1][mit->first].begin();
      std::advance(it, mit->second);
      if (it->idRange.first > up || it->idRange.second < lo) break;
      if (it->idRange.first > lo) toFillIn.push_back(EndPoint(row, lo, it->idRange.first - 1));
      record.push_back(it->idRange.second);
      lo = it->idRange.first;
      up = it->idRange.second;
      ++row;
    }
    //check if this row is fit for the query length
    if (row + 1 + k == query.length()) {
      while (lo <= up) {
        if (data[lo].length() == row + 1) {
          res.push_back(make_pair(lo, ed));
          if (res.size() >= topk) return true;
        } else break;
        ++lo;
      }
    }
    
    
    if (lo <= up) toFillIn.push_back(EndPoint(row, lo, up));

    // trace back record
    lo = record.back() + 1;
    record.pop_back();
    --row;
    while (!record.empty()) {
      if (record.back() >= lo) toFillIn.push_back(EndPoint(row, lo, record.back()));
      lo = record.back() + 1;
      record.pop_back();
      --row;
    }
    ++begin;
  }
    
  if (lo <= upp) toFillIn.push_back(EndPoint(row, lo, upp));
  return false;
}

// select next start position (diagonals[f], id, idx, shortest, lo). For example, when ed = 0 then id = 1. then this will return true.
bool NextStartPosition(const vector< vector<EndPoint> > & dia, int k, int idx[], int sh, int & stPos) {
  stPos = dia[k + sh][idx[sh + 1]].range.first;

  // i1 : next idx of shortest
  int i1 = (sh == 1 ? -1 : sh + 1);
  while (idx[i1 + 1] < dia[k + i1].size() && dia[k + i1][idx[i1 + 1]].range.second < stPos) idx[i1 + 1]++;
  if (idx[i1 + 1] >= dia[k + i1].size()) return false;
  else if (dia[k + i1][idx[i1 + 1]].range.first > stPos) return NextStartPosition(dia, k, idx, i1, stPos);

  int i2 = (i1 == 1 ? -1 : i1 + 1);
  while (idx[i2 + 1] < dia[k + i2].size() && dia[k + i2][idx[i2 + 1]].range.second < stPos) idx[i2 + 1]++;
  if (idx[i2 + 1] >= dia[k + i2].size()) return false;
  else if (dia[k + i2][idx[i2 + 1]].range.first > stPos) return NextStartPosition(dia, k, idx, i2, stPos);

  return true;
}





//return the k strings that satisfy for any string r in R (|R|=k), ED(r,q)<=ED(s,q) (s = set of string S - R)
void TopkSearch(const string & query, const unsigned topk) {
  int ed = -1;   //initialize the edit distance
  int idx[3];    //
  bool f = false;  //

  vector< vector<EndPoint> > diagonals[2];
    
    //int row, range<unsigned, unsigned>
  EndPoint edge(UNDEFROW, 0, data.size() - 1);
    
    //diagonal[f] has 3 endpoints each
  diagonals[f].assign(3, vector<EndPoint>(0));//initaliz the diagonal[f], make it has 3 vectors, and each vector has no EngPoint structure
  diagonals[f][0].push_back(edge);
  diagonals[f][1].push_back(EndPoint(-2, 0, data.size() - 1));
  diagonals[f][2].push_back(edge);

  
    
    
    //if the repository of R if smaller than k, then we have to continue the add edit distance each iteration and calculate.
   while(res.size() < topk) {
    ++ed;
    diagonals[!f].assign(2 * ed + 5, vector<EndPoint>(0)); //ini: diagonals[!f]:(5, no EndPoint right now);2nd time is 7, then 9,11...

    // notice that id and dia are a little tricky here.
    // the offsets of Diagonal[f] and Diagonal[!f] are different
    // offset of Diagonal[f] is dia + ed + 1
    // offset of Diagonal[!f] is dia + ed + 2
    // next using D[f][id-1], D[f][id] and D[f][id+1] to fill in D[!f][id+1]
      
       
       //the dia is the difference between the 2 strings such that it can be either negtive and positive. Obviously, if the ed is 0, the dia is definitly 0. if ed is 1, then the dia can be no bigger than 1 or -1.
       //the id is the map the dia to array index
       for (int dia = -1 * ed; dia <= ed; dia++) {
           //initialize the id, shortest and idx[] for this dia in this edit distance
      int id = dia + ed + 1;// when ed = 0, dia = 0, id=1; |**| ed = 1, dia = -1(ini), id = 1 | dia = 0, id = 2 | dia = 1, id = 3; |**| ed=2, dia=-2, id=1 | dia = -1,id=2|dia=0,id=3..|dia=2,id=5
      
      int lo, up, row, shortest = 0; //lo means the lower bound of range, up = upper bound, row means how deep the node, shortest means the length range closest to query
      idx[0] = idx[1] = idx[2] = 0;//store the number of endpoint for each length of extenstions
       
      // test if there is 'next' string id range to search
      // shortest is a value of -1, 0, 1; shortest + 1 maps to 0,1,2
      // id + shortest is map to id - 1, id, id + 1 (actually, no smaller than 0)
      while (idx[shortest + 1] < diagonals[f][id + shortest].size()  //idx[1] < dia[f][1+0].size() && NeStPo(dia[f], 1, 0, 0, 0)
          && NextStartPosition(diagonals[f], id, idx, shortest, lo))  {

        // select a maximum row
        // row = maximum(diagonals[f][id - 1][idx[0]].row,
        // diagonals[f][id][idx[1]].row + 1, diagonals[f][id + 1][idx[2]].row + 1);
        row = diagonals[f][id - 1][idx[0]].row <= diagonals[f][id][idx[1]].row + 1
            ? (diagonals[f][id + 1][idx[2]].row <= diagonals[f][id][idx[1]].row
            ? diagonals[f][id][idx[1]].row + 1 : diagonals[f][id + 1][idx[2]].row + 1)//divide the 2 options here
            : (diagonals[f][id - 1][idx[0]].row <= diagonals[f][id + 1][idx[2]].row + 1
            ? diagonals[f][id + 1][idx[2]].row + 1 : diagonals[f][id - 1][idx[0]].row);

        // select a minimum end of range
        shortest = (diagonals[f][id - 1][idx[0]].range.second <= diagonals[f][id][idx[1]].range.second 
            ? (diagonals[f][id - 1][idx[0]].range.second <= diagonals[f][id + 1][idx[2]].range.second ? -1 : 1)
            : (diagonals[f][id][idx[1]].range.second <= diagonals[f][id + 1][idx[2]].range.second ? 0 : 1));
        //find the upper bound using updated shortest
        up = diagonals[f][id + shortest][idx[shortest + 1]].range.second;
          
        //add one endpoint of max row and min range
        idx[shortest + 1]++;
        //basic computation is done here, below is to compare then modify return the corresponding string
          
        if (row + 1 + dia > query.length()) continue;//if the deepth of row+difference is beyond the query length, then goto next dia
        //or if they are equal
        else if (row + 1 + dia == query.length()) {
          for (unsigned x = lo; x <= up; x++) {
            if (data[x].length() == row + 1) {
              res.push_back(make_pair(x, ed));
                if (res.size() >= topk) return; //check if this reaches the k
            }
          }
            
          //if row + 1 + dia == query.length(). put it into diagonals[!f]. this will be used as input of next iteration
          diagonals[!f][id + 1].push_back(EndPoint(row, lo, up));
            
            //the rest of it will handle strings whose length is not match row+1
            //if row + 1 + dia < query.length() and
        } else if (up - lo + 1 < 4) {
             //to store present value of row
            int curRow = row;
            //first, test the lo bound
            //if the char at row+1 position of data[lo] is the same with query[row+dia+1], just add row. go to a deeper node
            //and if this addition makes the length match, then we found an answer and add it to result set
            while(row + 1 + dia < query.length() && row + 1 < data[lo].length() && data[lo][row + 1] == query[row + dia + 1]) row++;
          if (row + 1 + dia == query.length() && row + 1 == data[lo].length()) {
            res.push_back(make_pair(lo, ed));
            if (res.size() >= topk) return;//return if it reaches
          }
          //if not, put it as endpoint
          else diagonals[!f][id + 1].push_back(EndPoint(row, lo, lo));
            
            //if the length of query and the length of data string with lower bound are not equal to row +1 +dia and row + 1 respectively
            //notice that prefixLen occurred in main function, it recorded the location of different char in each string compared to last string
          for (unsigned x = lo + 1; x <= up; x++) {//change upper range of this endpoint to x
            if (prefixLen[x] > row + 1) diagonals[!f][id + 1].back().range.second = x;
            else if (prefixLen[x] <= curRow) {
              row = curRow;
              while (row + 1 + dia < query.length() && row + 1 < data[x].length() && data[x][row + 1] == query[row + dia + 1])
                row++;
              if (row + 1 + dia == query.length() && row + 1 == data[x].length()) {
                res.push_back(make_pair(x, ed));
                if (res.size() >= topk) return;//return if it reaches
              } else diagonals[!f][id + 1].push_back(EndPoint(row, x, x));
              //if the prefixlen is smaller than row and larger than curRow
            } else if (prefixLen[x] <= row && prefixLen[x] > curRow) {
              row = prefixLen[x] - 1;
              diagonals[!f][id + 1].push_back(EndPoint(row, x, x));
            } else {
              while (row + 1 + dia < query.length() && row + 1 < data[x].length() && data[x][row + 1] == query[row + dia + 1])
                row++;
              if (row + 1 + dia == query.length() && row + 1 == data[x].length()) {
                res.push_back(make_pair(x, ed));
                if (res.size() >= topk) return; // if it reaches
              } else if (diagonals[!f][id+1].back().row == row) diagonals[!f][id + 1].back().range.second = x;
              else diagonals[!f][id + 1].push_back(EndPoint(row, x, x));
            }
          }
        }
        //if row + 1 + dia < query.length() and up - lo + 1 >= 4
        else if (RollDown(diagonals[!f][id + 1], dia, row, lo, up, query, ed, topk)) return;
      }
    }

       
    //no more pivotal entries in this edit distance, using the !f. the f and !f
    diagonals[!f][0].push_back(edge);
    diagonals[!f][1].push_back(edge);
    diagonals[!f][2 * ed + 3].push_back(edge);
    diagonals[!f][2 * ed + 4].push_back(edge);
    f = !f;
  }
}
