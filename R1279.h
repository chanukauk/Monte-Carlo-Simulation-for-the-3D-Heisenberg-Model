#ifndef PSIMAG_R1279_H
#define PSIMAG_R1279_H

/* Shift register random number generator */

/*
   Encapsulated and made compliant to psimag
   interface by : H. K. Lee, 19 June 2002 
*/


/* original c version was converted from a fortran */
/* version by WooSeap Kwak                         */


#include <iostream>
#include <stdexcept>

class R1279 {

  public:
    typedef double result_type;
    typedef int  seed_type;
    static const bool has_fixed_range = true; 
    static result_type min() { return 0; }
    static result_type max() { return 1; }

  private:
    static const int mask = 1;
    static const int nab3 = 101280;
    unsigned int     num;
    int              ranvec[nab3];
    result_type      ranvec2[nab3];

  public:
    //--------------------------------------------------------------
    explicit R1279(seed_type iseed=1) {
      seed(iseed);
    }

    //--------------------------------------------------------------
    void seed(seed_type iseed) { 
      if(iseed<=0) iseed = -iseed+1;
      num = nab3 - 1280;
      ranini(iseed);
    }
    //--------------------------------------------------------------
    result_type operator()() {
      result_type temp;
      temp = randa(num);
      num++;
      return temp;
    }
    //--------------------------------------------------------------
    result_type operator()(result_type n) {
      return (*this)()*n;
    }
    //--------------------------------------------------------------
    friend std::istream& operator >> (std::istream& is, R1279& rnd) {
      seed_type iseed;
      is >> iseed; rnd.seed(iseed); return is;
    }
    //--------------------------------------------------------------
    friend bool operator == (const R1279& a, const R1279& b) {

      if(a.num!=b.num) return false; 

      for(unsigned int i=0;i<nab3;i++) {
	if(a.ranvec[i]!=b.ranvec[i]) return false; 
      }
      return true;
    }
    //--------------------------------------------------------------

  private:
    //--------------------------------------------------------------
    result_type epsilon() { return .0000000004656612875; }
    //--------------------------------------------------------------
    void ranini(int iseed) {

      int       imod;
      int       imax;
      result_type      rmod, pmod, dmaxi;
      int       tmask, jt;

      imax  = 2147483647;
      dmaxi = 1.0/2147483647.0;
      rmod  = static_cast<result_type>(iseed);
      pmod  = static_cast<result_type>(imax);


      for(int i=0;i<1000;i++) {
	rmod = rmod * 16807.0;
	imod = (int)(rmod * dmaxi);
	rmod = rmod - pmod * imod;
      }

      for(int i=0; i<1279;i++) {
	ranvec[i]=0;
	for(int j=0;j<31;j++) {
	  for(int k=0;k<36;k++) {
	    rmod = rmod * 16807.0;
	    imod = (int)(rmod * dmaxi);
	    rmod = rmod - pmod * imod;
	  }
	  rmod = rmod * 16807.0;
	  imod = (int)(rmod * dmaxi);
	  rmod = rmod - pmod * imod;
	  if( rmod > ( 0.5 * pmod ) ) {
	    tmask = mask;
	    jt=0;
	    while(jt<j) {	
	      tmask <<=1;
	      jt++;
	    }
	    ranvec[i]|=tmask;
	  }
	}
      }
      random2(1000);
    }
    //--------------------------------------------------------------
    void random2(unsigned int number) {

      for(unsigned int i=0;i<number;i+=4) {
	ranvec[i+1279]=ranvec[i  ]^ranvec[i+216];
	ranvec[i+1280]=ranvec[i+1]^ranvec[i+217];
	ranvec[i+1281]=ranvec[i+2]^ranvec[i+218];
	ranvec[i+1282]=ranvec[i+3]^ranvec[i+219];
      }
      for(unsigned int i=0;i<1276;i+=4) {
	ranvec[i  ] = ranvec[i  +number];
	ranvec[i+1] = ranvec[i+1+number];
	ranvec[i+2] = ranvec[i+2+number];
	ranvec[i+3] = ranvec[i+3+number];
      }
      ranvec[1276]=ranvec[1276+number];
      ranvec[1277]=ranvec[1277+number];
      ranvec[1278]=ranvec[1278+number];
      for(unsigned int i=0;i<number;i++) {
	ranvec2[i]=static_cast<result_type>(ranvec[i]*epsilon());
      }
    }
    //--------------------------------------------------------------
    result_type randa(unsigned int number) {
      result_type temp;
      if(number >= (nab3-1280) ) {
	random2(nab3-1280);
	temp = ranvec2[0];
	num = 2;
      }
      else {
	temp = ranvec2[number-1];
	num++;
      }
      return temp;
    }
};



#endif


