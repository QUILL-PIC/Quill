#include <R.h>
#include <Rinternals.h>
#include <iostream>
#include <cstring>

using namespace std;

// example of R function written in C
extern "C"
SEXP foo(SEXP i) {
    SEXP res = PROTECT(allocVector(INTSXP, 1));
    INTEGER(res)[0] = 2 * INTEGER(i)[0] + 1;
    UNPROTECT(1);
    return res;
}

// chameleon-core (libchameleon) functions
extern "C" void configure(char* charlog, bool flag = false);
extern "C" double get(char* charname);
extern "C" double* read2d(char* charfilename, char* charplane);

// bindings for R; "fr" stands for "from R"

const int max_str_length = 100;

void sexpcpy(char* cstr, SEXP string) {
    int l = strlen(CHAR(STRING_ELT(string, 0)));
    if (l < max_str_length) {
        strcpy(cstr, CHAR(STRING_ELT(string, 0)));
    } else {
        cerr << "libchameleon-R: sexpcpy: ERROR: length of input string is too"
            "high, void string is returned." << endl;
        char c[] = "";
        strcpy(cstr, c);
    }
}

SEXP zero() {
    SEXP res = PROTECT(allocVector(INTSXP, 1));
    INTEGER(res)[0] = 0;
    UNPROTECT(1);
    return res;
}

extern "C"
SEXP configure_fr(SEXP log, SEXP flag) {
    char charlog[max_str_length];
    sexpcpy(charlog, log);
    configure(charlog, LOGICAL(flag)[0]);
    //
    return zero();
}

extern "C"
SEXP get_fr(SEXP name) {
    char charname[max_str_length];
    sexpcpy(charname, name);
    SEXP value = PROTECT(allocVector(REALSXP, 1));
    REAL(value)[0] = get(charname);
    UNPROTECT(1);
    return value;
}

extern "C"
SEXP read_fr(SEXP filename, SEXP plane) {
    char charfilename[max_str_length];
    sexpcpy(charfilename, filename);
    char charplane[max_str_length];
    sexpcpy(charplane, plane);
    double* a = read2d(charfilename, charplane);
    long n;
    string strplane(charplane);
    if (strplane == "xy") {
        char nx[] = "nx";
        char ny[] = "ny";
        n = get(nx) * get(ny);
    } else if (strplane == "xz") {
        char nx[] = "nx";
        char nz[] = "nz";
        n = get(nx) * get(nz);
    } else {
        char ny[] = "ny";
        char nz[] = "nz";
        n = get(ny) * get(nz);
    }
    SEXP b = PROTECT(allocVector(REALSXP, n));
    for (long i = 0; i < n; ++i)
        REAL(b)[i] = a[i];
    UNPROTECT(1);
    delete[] a;
    return b;
}
