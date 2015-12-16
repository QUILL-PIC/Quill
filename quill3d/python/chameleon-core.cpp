#include <fstream>
#include <iostream>
#include <map>

using namespace std;

static map<string, double> conf;

// fills up global variable conf with values from log file,
// prints variable name and value if flag is true
extern "C"
void configure(const char* charlog, bool flag = false) {
    string log(charlog); // copies the null-terminated character sequence
    fstream fs(log, ios_base::in);
    if (fs) {
        string name;
        double value;
        while (fs >> name && name[0] != '#') { /* # is the terminating symbol in
                                                  log */
            if (fs >> value) {
                conf[name] = value;
                if (flag)
                    cout << name << " = " << value << '\n';
            } else {
                // this can happen if value is a string, not a double; then
                // error state flags are cleared and new attempt of name reading
                // occurs
                fs.clear();
            }
        }
    } else {
        cerr << "libchameleon: configure: ERROR: no such file, " << log <<
            endl;
    }
}

// returns value from conf that corresponds to charname
extern "C"
double get(const char* charname) {
    string name(charname);
    return conf[name];
}

// reads 2D data
extern "C"
double* read2d(const char* charfilename, const char* charplane) {
    string filename(charfilename);
    string plane(charplane);
    long nx = conf["nx"];
    long ny = conf["ny"];
    long nz = conf["nz"];
    long n = nx * (ny + nz) + ny * nz;
    double* a;
    fstream fs(filename, ios_base::in);
    if (fs) {
        if (conf["output_mode"] == 0) {
            // text mode
            a = new double[n];
            long i = 0;
            while (fs >> a[i] && i < n)
                ++i;
        } else if (conf["output_mode"] == 1) {
            // binary mode
            fs.close();
            fs.open(filename, ios_base::in | ios_base::binary);
            a = new double[n];
            fs.read(reinterpret_cast<char*>(a), sizeof(double) * n);
        } else {
            cerr << "libchameleon: read:: ERROR: unknown output_mode" << endl;
            return nullptr;
        }
    } else {
        cerr << "libchameleon: read: ERROR: no such file, " << filename <<
            endl;
        return nullptr;
    }
    fs.close();
    double* b;
    if (plane == "xy") {
        b = new double[nx * ny];
        for (long i = 0; i < nx; ++i)
            for (long j = 0; j < ny; ++j)
                b[i * ny + j] = a[i * (ny + nz) + j];
    } else if (plane == "xz") {
        b = new double[nx * nz];
        for (long i = 0; i < nx; ++i)
            for (long j = 0; j < nz; ++j)
                b[i * nz + j] = a[i * (ny + nz) + ny + j];
    } else if (plane == "yz") {
        b = new double[ny * nz];
        for (long i = 0; i < ny; ++i)
            for (long j = 0; j < nz; ++j)
                b[i * nz + j] = a[nx * (ny + nz) + i * nz + j];
    } else {
        cerr << "libchameleon: read: possible plane values are xy, xz and yz,"
            "given: " << plane << endl;
        delete[] a;
        return nullptr;
    }
    delete[] a;
    return b;
}
