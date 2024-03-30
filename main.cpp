// --------------------------------------------------------------------------
// -----------------------------Include--------------------------------------
// --------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <float.h>
#include <cmath>
#include <iostream>
#include <lapack.h>
#include <cholmod.h>
#include <SuiteSparse_config.h>
#include <SuiteSparseQR_C.h>

// --------------------------------------------------------------------------
// ---------------------------------宏---------------------------------------
// --------------------------------------------------------------------------
// 宏定义 Linux
#define LINUX

#ifdef LINUX
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#elif defined(WINDOWS)
#include <windows.h>
#endif

// ---------------------------------------------------------------------------
// -----------------------------Namespace-------------------------------------
// ---------------------------------------------------------------------------
//using namespace Eigen;
using namespace std;

// ---------------------------------------------------------------------------
// -----------------------------结构体定义--------------------------------------
// ---------------------------------------------------------------------------
// 假设 SitesInfo 结构体
typedef struct {
    char** name;
    int* doy;
    double** coor;
    double* RDCB_REF;
} SitesInfo;

typedef struct {
    //P1，P2，L1，L2均为二维double数组
    double **P1;
    double **P2;
    double **L1;
    double **L2;
} Obs;

typedef struct {
    double** value;
    int* doy;
} SDCB_REF;

// -------------------------------------------------------------------------
// -----------------------------全局变量--------------------------------------
// -------------------------------------------------------------------------
//const char* r_ipath = "/home/jason/projects/M_DCB_C/RINEX_files";// rinex
const char* r_ipath = "/home/jason/projects/M_DCB_C/data/2021_029_data/obs";// rinex
//const char* r_opath = "/home/jason/projects/M_DCB_C/RINEX_output_files";//sprintf(sav_filename, "D:\\projects\\M_DCB\\RINEX_output_files\\observation_%s.csv", psitesInfo->name[i]);
const char* r_opath = "/home/jason/projects/M_DCB_C/data/2021_029_data/obs_output";
const char* s_ipath = "/home/jason/projects/M_DCB_C/SP3_files";// sp3
const char* s_opath = "/home/jason/projects/M_DCB_C/SP3_output_files";//
const char* i_ipath = "/home/jason/projects/M_DCB_C/IONEX_files";// ionex
const char* i_opath = "/home/jason/projects/M_DCB_C/IONEX_output_files";// ionex
const char* m_p4_path = "/home/jason/projects/M_DCB_C/M_P4";// P4
const char* m_result_path = "/home/jason/projects/M_DCB_C/M_Result";// result
int lim = 10;// el
int order = 4;// order
int r_file_num;
double* DCB_rec;
double f1 = 1575.42e6;
double f2 = 1227.6e6;
double pi = 3.14159265358979323846;
double c = 299792458;
int period = 2;
SitesInfo sitesInfo;
Obs obs;
SDCB_REF sDCB_REF;

const char* OBS_type_1 = "C2I";
const char* OBS_type_2 = "C6I";
const char* Constellation = "C";
int N_SAT = 45;

// --------------------------------------------------------------------------
// -----------------------------函数声明--------------------------------------
// --------------------------------------------------------------------------
int countFilesInDirectory(const char *folderPath);
void read_rinex(const char* r_ipath, const char* r_opath, SitesInfo *psitesInfo, Obs *pobs);
void writeArrayToFileWithLabels(FILE *file, double **array, const char* dataType);
void read_sp3(const char* s_ipath, const char* s_opath);
void parse_sp3(char* sp3_file, double ***xyz);
void interplotation(double ***pre_xyz, double ***cur_xyz, double ***next_xyz, double ***sate_xyz);
void extend_matrix(int dimension, double ***pre_xyz, double ***cur_xyz, double ***next_xyz, double ***xyz_etd);
double* interp_lag(double* x, double* y, double* x0);
int* GWeek_2_DOY(int G_Week, int Day_of_Week);
void saveSp3ToCSV(double*** state_xyz, char* filename);
void read_ionex(const char* i_ipath, const char* i_opath, SitesInfo *psitesInfo, SDCB_REF sdcb_ref);
int find_ionex_index(int doy, char** ionex_filenames, int ionex_file_num);
int find_rinex_indexs(int* index, int doy);
void read_single_ionex(const char* ionex_file, int i, char** sites_name, int doy_num, SDCB_REF sDCB_REF, double* DCB_rec);
void get_smoothed_P4(SitesInfo sitesInfo, double z_threshold, int flag);
void writeObsToFile(const char* filename, Obs* obs, int p1Rows, int p2Rows, int l1Rows, int l2Rows, int cols);
void readObsFromFile(const char* filename, Obs* obs);
void save_sp3_To_Bin_File(double ***sate_xyz, int dim1, int dim2, int dim3, const char *filename);
void load_Sp3_From_Bin_File(double ****sate_xyz, int dim1, int dim2, int dim3, const char *filename);
double* XYZ2BLH(double x, double y, double z);
double* get_EA(double sx, double sy, double sz, double x, double y, double z);
void cutobs(Obs *obs_temp, double*** sate_xyz, double sx, double sy, double sz, double lim);
double** pre_pro(Obs* obs_temp);
int** Get_arc(double** L6, int prn,int* parc_n);
void deleteRowAndReplace(int ***arc, int *n, int row);
double mean(double** array, int start, int end, int prn);
void save_P4_To_Bin_File(double **P4, int dim1, int dim2, const char *filename);
void load_P4_From_Bin_File(double ***P4, int dim1, int dim2, const char *filename);
char** getFileNamesInDirectory(const char *path, char **fileNames, int num);
void DCB_Estimation();
int** ConvertVectorTo2DMatrix(int *arr, int size, int *returnSize);
void InsertRowInArc(int ***arc, int *n, int row, int k, int k_shift);
void removeSameElementAndReturn_Char(char** vector, char*** new_vector, int vector_size, int* new_size);
void removeSameElementAndReturn_Int(int* vector, int** new_vector, int vector_size, int* new_size);
char* get_load_sp3_pathname(int doy);
char** get_doy_P4_filepath(int doy, int* n_r);
void get_DCB(double** DCB_R, int* DCB_R_size, double** DCB_S, int* DCB_S_size, double** IONC, int* IONC_size, int doy, double*** sate_xyz, SitesInfo sitesInfo, SDCB_REF sDCB_REF, int order);
double* get_IPP(double E, double A, double sb, double sl, double IPPz, double t_r);
double* legendre(int n, double x);
double norm(int n, int m);
double factorial(int n);
void get_Coef(double** M_col, int st, int ed, double b, double s, int order);
void get_Matrix(double** P4, double** x, double** y, double**z, double sx, double sy, double sz, int n_s, int n_r, int ith, int order, double*** sN, double** sL, int est_num, int* M_row_number, int* l_size_number);
double* get_legendre(int n, double x);
void addRowTo2DMatrixNewLine(double ***M, double **M_sol, int *numRows, int numCols);
void appendDoubleToVector(double **arr, int *size, double value);
void addMartixDownToAnother(double ***B, double*** sN, int est_num, int* B_row_number, int* l_row_number, int sN_row_number);
void addVectorToVector(double** l, double** sL, int* l_row_number, int sL_size);
double* getDCB_R(double* R, int n_r);
double* getDCB_S(double* R, int n_r, int n_s);
double* getIONC(double* R, int n_r, int n_s);
void write_vector_to_file(char *filename, double* vector, int size);
char* generate_final_result_file_pathname(int doy, const char* type, const char* m_result_path);
void malloc_Double_2D(double*** array, int row, int col);
void malloc_Double_Vector(double** vector, int size);
void malloc_Char_Vector(char*** vector, int size, int string_size);
void malloc_Int_Vector(int** vector, int size);
void malloc_double_3D(double**** array, int dim1, int dim2, int dim3);
void malloc_Int_2D(int*** array, int row, int col);
void sort_filenames(char *files[], int count);
int compare(const void *a, const void *b);
int strnatcmp(const char *a, const char *b);
double* solve_least_squares(double** B, double* l, int B_row_num, int B_col_num, int l_row_num) ;
void read_rinex_v304(const char* r_ipath, const char* r_opath, SitesInfo *psitesInfo, Obs *pobs,const char* OBS_type_1, const char* OBS_type_2, const char* Constellation);


// 不同平台的函数声明
#ifdef LINUX
void CreateFolderIfNotExist_Linux(const char* folderPath);
#elif defined(WINDOWS)
void CreateFolderIfNotExist_WINDOWS(const char* folderPath);
#endif

// -----------------------------------------------------------------------
// -----------------------------Main--------------------------------------
// -----------------------------------------------------------------------
int main() {
    setbuf(stdout, 0); //开启debug时console也有输出信息
    // 打印基本信息
    printf("elevation angle threshold(unit:degree)(10 degree is recommended): %d\n", lim);
    printf("the order of spheric harmonic function (4 order is recommended): %d\n", order);
    printf("MDCB(multi-stations) starts running!\n");

    // 创建文件夹，linux和windows创建方法不一样
#ifdef LINUX
    CreateFolderIfNotExist_Linux(r_opath);
    CreateFolderIfNotExist_Linux(s_opath);
    CreateFolderIfNotExist_Linux(i_opath);
    CreateFolderIfNotExist_Linux(m_p4_path);
    CreateFolderIfNotExist_Linux(m_result_path);
#elif defined(WINDOWS)
    CreateFolderIfNotExist_WINDOWS(r_opath);
    CreateFolderIfNotExist_WINDOWS(s_opath);
    CreateFolderIfNotExist_WINDOWS(i_opath);
    CreateFolderIfNotExist_WINDOWS(m_p4_path);
    CreateFolderIfNotExist_WINDOWS(m_result_path);
#endif

    // 为obs结构体内的数组分配内存，每个变量均需要一个(2880,N_SAT)尺寸的二维double数组
    malloc_Double_2D(&obs.P1, 2880, N_SAT);  // 分配内存空间给 P1
    malloc_Double_2D(&obs.P2, 2880, N_SAT); // 分配内存空间给 P2
    malloc_Double_2D(&obs.L1, 2880, N_SAT); // 分配内存空间给 L1
    malloc_Double_2D(&obs.L2, 2880, N_SAT); // 分配内存空间给 L2

    // Step one----read rinex files
    read_rinex_v304(r_ipath, r_opath, &sitesInfo, &obs, OBS_type_1, OBS_type_2, Constellation);
//    read_rinex(r_ipath, r_opath, &sitesInfo, &obs);

    // Step two----Read SP3 files
    read_sp3(s_ipath, s_opath);

    // Step three----Read ionex files
    //read_ionex(i_ipath, i_opath, &sitesInfo, sDCB_REF);

    // Step four----Get smoothed P4
    get_smoothed_P4(sitesInfo, (lim*pi/180), 0);

    // Step Five----DCB Estimation
    DCB_Estimation();

    return 0;
}

// ------------------------------------------------------------------------
// -----------------------------函数区--------------------------------------
// ------------------------------------------------------------------------

// 功能描述：读取rinex文件，把P1,P2,L1,L2观测数据存为二进制文件，测站相关信息存入psitesInfo结构体中
// 入参：r_ipath —— 存放rinex文件夹路径； r_opath —— 输出观测数据文件的文件夹路径；
//      psitesInfo —— 事先定义好的存放测站信息的结构体； pobs —— 事先定义好的存放观测信息的结构体
// 出参：无
void read_rinex(const char* r_ipath, const char* r_opath, SitesInfo *psitesInfo, Obs *pobs) {
    printf("----------Step 1: starting read rinex files!----------\n");
    DIR *dir;
    FILE *file;
    struct dirent *entry;

    // 打开文件夹
    dir = opendir(r_ipath);
    if (dir == NULL) {
        perror("无法打开文件夹");
        exit(EXIT_FAILURE);
    }

    printf("Starting Reading Rinex Files!\n");

    r_file_num = countFilesInDirectory(r_ipath);    // 计算rinex文件数
    char** fileNames;   //存放rinex文件夹下的所有文件名
    malloc_Char_Vector(&fileNames, r_file_num, 256);    //为fileNames分配空间
    getFileNamesInDirectory(r_ipath, fileNames, r_file_num);    //获取rinex文件夹下的所有文件名
    sort_filenames(fileNames, r_file_num);  //对fileNames文件名进行自然排序

    // 为 SitesInfo 结构体内的变量分配内存
    malloc_Char_Vector(&psitesInfo->name, r_file_num, 256);
    malloc_Int_Vector(&psitesInfo->doy, r_file_num);
    malloc_Double_2D(&psitesInfo->coor, r_file_num, 3);
    malloc_Double_Vector(&psitesInfo->RDCB_REF, r_file_num);

    // 为obs结构体内的数组分配内存，每个变量均需要一个(2880,32)尺寸的二维double数组
    malloc_Double_2D(&pobs->P1, 2880, N_SAT);  // 分配内存空间给 P1
    malloc_Double_2D(&pobs->P2, 2880, N_SAT); // 分配内存空间给 P2
    malloc_Double_2D(&pobs->L1, 2880, N_SAT); // 分配内存空间给 L1
    malloc_Double_2D(&pobs->L2, 2880, N_SAT); // 分配内存空间给 L2

    char doy[4];    //用于存放文件名中的doy部分
    char sub_str[3];    //用于存放文件名中的
    // 循环rinex文件名，取出测站名和doy存入psitesInfo里（读sitesInfo）
    for (int i = 0; i < r_file_num; i++) {
        strcpy(psitesInfo->name[i], fileNames[i]);  //fileNames不处理，直接存入name
        //将fileNames[i]的第10到11个字符赋值给sub_str，再将fileNames[i]的第5到7个字符赋值给doy
        //例如，文件名叫gope0010.10o，sub_str就是10，doy就是001
        strncpy(sub_str, &fileNames[i][9], 2);
        strncpy(doy, &fileNames[i][4], 3);
        sub_str[2] = '\0';  //这一步必须加入字符串结束标志，否则会出问题
        doy[3] = '\0';  //这一步必须加入字符串结束标志，否则会出问题
        psitesInfo->doy[i] = 1000 * atoi(sub_str) + atoi(doy);  //计算真正的doy，存入doy
    }
    // 释放fileNames的内存
    //free(fileNames);
//    //打印psitesInfo中存好的测站名和doy，用于debug
//    for (int i = 0; i < r_file_num; i++) {
//        printf("File No.%d :%s, DOY: %d\n", i + 1, psitesInfo->name[i], psitesInfo->doy[i]);
//    }

    char line[256]; // 假设读文件过程中一行最多包含256个字符
    char filename[256]; //待读的文件名称
    int obst_n; char **obst;
    double h; double m; double s; int ep = 0;
    int *loc;

    // 循环读取每个rinex（读obs观测值）
    for (int i = 0; i < r_file_num; i++) {
        //printf("Reading File No.%d :%s\n", i + 1, psitesInfo->name[i]);

        //拼接部分读取文件的路径
        strcpy(filename, r_ipath);
#ifdef LINUX
        strcat(filename, "/");  //Linux文件系统和windows文件路径符号不一样
#elif defined(WINDOWS)
        strcat(filename, "\\");
#endif
        strcat(filename, psitesInfo->name[i]);  //拼接待读文件完整路径
        file = fopen(filename, "r");    // 打开文件
        if (file == NULL) {
            perror("Can't Open File");
        }

        // 逐行读取
        while (fgets(line, sizeof(line), file)) {
            // -----get receivers coordinates-----------
            if (strstr(line, "APPROX POSITION XYZ") != NULL) {
                double x, y, z;
                // 将行解析成三个双精度浮点数
                if (sscanf(line, "%lf %lf %lf", &x, &y, &z) == 3) {
                    psitesInfo->coor[i][0] = x;
                    psitesInfo->coor[i][1] = y;
                    psitesInfo->coor[i][2] = z;
                }
            }
            //-----get the GPS observables' types------
            if (strstr(line, "# / TYPES OF OBSERV") != NULL){
                if(sscanf(line, "%d", &obst_n) == 1){
                    if (obst_n > 9){
                        printf("obst_n > 9!!!!!");
                    }
                    else{
                        malloc_Char_Vector(&obst, obst_n, 4);
                        int offset = 10;
                        for (int j = 0; j < obst_n; j++){
                            if (sscanf(line + offset, "%2s", obst[j]) == 1) {
                                // 成功解析一个观测名称，更新偏移量，准备解析下一个观测名称
                                offset += 6; // 每个观测名称占据4个字符，包括空格
                            } else {
                                // 解析失败，处理异常情况
                                printf("Failed to parse observation type at index %d\n", j);
                                break;
                            }
                        }

                        malloc_Int_Vector(&loc, obst_n);
                        for (int j = 0; j < obst_n; j++){
                            //比较obst[j]和"L1"字符串是否相同，如果相同，则给loc[0]赋值为j
                            if (strcmp(obst[j], "L1") == 0){
                                loc[0] = j;
                            }else if (strcmp(obst[j], "L2") == 0) {
                                loc[1] = j;
                            }else if (strcmp(obst[j], "P1") == 0) {
                                loc[2] = j;
                            }else if (strcmp(obst[j], "P2") == 0) {
                                loc[3] = j;
                            }
                        }
                        if (loc[2] == 0 || loc[3] == 0){
                            printf("P1 and P2 not found!\n");
                        }
                    }
                    // 释放obst和loc的内存
                    //free(obst);
                }
            }
            //----start get GPS observables------------
            if (strstr(line, "END OF HEADER") != NULL) {
                while (1){
                    //读取下一行
                    fgets(line, sizeof(line), file);
                    // 检查是否到达文件结尾
                    if (feof(file)) {
                        // 到达文件结尾，跳出循环
                        break;
                    }
                    //如果line长度大于32个个字符，并且line的第33个字符是"G"或"R"或空格，则执行下面的语句
                    if (strlen(line) > 32 && (line[32] == 'G' || line[32] == 'R')){
                        //将line的第11到第12个字符转换成double类型，再赋值给h
                        h = 0; m = 0; s = 0;
                        if (sscanf(line + 10, "%2lf", &h) == 1){
                        }
                        if (sscanf(line + 13, "%2lf", &m) == 1){
                        }
                        if (sscanf(line + 16, "%10lf", &s) == 1){
                        }
                        ep = h*120 + m*2 + s/30 + 1;
                    }else{
                        continue;
                    }
                    //---get satellite number----

                    int nsat; int *sv_G;
                    //取出line的第31到32个字符，转换成int类型，再赋值给nsat
                    if (sscanf(line + 30, "%2d", &nsat) == 1){
                    }
                    //给sv_G申请nsat个int类型的内存，并将数组中每个值赋值为0
                    malloc_Int_Vector(&sv_G, nsat);
                    if (nsat>12){
                        for(int j = 0; j < 12; j++){
                            if (line[29+3*(j+1)] == 'G'){
                                //取出line的第31+3*j到32+3*j个字符，转换成int类型，再赋值给sv_G[j]
                                if (sscanf(line + 30 + 3*(j+1), "%2d", &sv_G[j]) == 1){
                                }
                            }
                        }
                        //读取下一行
                        fgets(line, sizeof(line), file);
                        if (nsat<25){
                            for (int j = 0; j < nsat-12; j++){
                                if (line[29+3*(j+1)] == 'G' || line[32] == ' '){
                                    if (sscanf(line + 30 + 3*(j+1), "%2d", &sv_G[j+12]) == 1){
                                    }
                                }
                            }
                        }
                    }else{
                        for(int j = 0; j < nsat; j++){
                            if (line[29+3*(j+1)] == 'G' || line[32] == ' '){
                                //取出line的第32+3*j到32+3*j+2个字符，转换成int类型，再赋值给sv_G[j]
                                if (sscanf(line + 30 + 3*(j+1), "%2d", &sv_G[j]) == 1){
                                }
                            }
                        }
                    }

                    //---get the observations----
                    for (int j = 0; j < nsat; j++){
                        fgets(line, sizeof(line), file);

                        double *obs_temp;
                        malloc_Double_Vector(&obs_temp, obst_n);

                        if (obst_n>5){
                            if (sv_G[j] == 0){
                                fgets(line, sizeof(line), file);
                                continue;
                            }
                            for (int k = 0; k < 5; k++){//
                                //判断line的长度，如果大于16*(k+1)-3，则执行下面的语句
                                if (strlen(line) > 16*(k+1)-3){
                                    //strcmp(line(16*j-15:16*j-2),'              '),如果line的第15*(k+1)-15到第15*(k+1)-2个字符是空格，则执行下面的语句
                                    if (strncmp(line + 16*k, "              ", 14) == 0){
                                        continue;
                                    }
                                    //取出line的第15*(j+1)-15至第15*(j+1)-2个字符，转换成double类型，再赋值给obs_temp[j]
                                    if (sscanf(line + 16*k, "%14lf", &obs_temp[k]) == 1){
                                    }
                                }
                            }
                            fgets(line, sizeof(line), file);
                            for (int k = 0; k < obst_n-5; k++) {
                                // 判断line的长度，如果大于16*(k+1)-3，则执行下面的语句
                                if (strlen(line) > 16*(k+1)-3) {
                                    // strcmp(line + 16*k, "              "), 如果line的第16*k至第16*k+14个字符是空格，则执行下面的语句
                                    if (strcmp(line + 16*k, "              ") == 0) {
                                        continue;
                                    }
                                    // 取出line的第16*k至第16*k+14个字符，转换成double类型，再赋值给obs_temp[k+5]
                                    // 注意：C语言中数组索引从0开始，因此对应MATLAB中的j+5需要在C语言中调整为k+5
                                    if (sscanf(line + 16*k, "%14lf", &obs_temp[k+5]) == 1) {
                                        // 成功转换字符串为double并赋值
                                    }
                                }
                            }
                        }else{
                            if (sv_G[j] == 0){
                                fgets(line, sizeof(line), file);
                                continue;
                            }
                            for (int k = 0; k < obst_n; k++){
                                //判断line的长度，如果大于16*(k+1)-3，则执行下面的语句
                                if (strlen(line) > 16*(k+1)-3){
                                    //strcmp(line(16*j-15:16*j-2),'              '),如果line的第15*(k+1)-15到第15*(k+1)-2个字符是空格，则执行下面的语句
                                    if (strncmp(line + 16*k, "              ",14) == 0){
                                        continue;
                                    }
                                    //取出line的第15*(j+1)-15至第15*(j+1)-2个字符，转换成double类型，再赋值给obs_temp[j]
                                    if (sscanf(line + 16*k, "%14lf", &obs_temp[k]) == 1){
                                    }
                                }
                            }
                        }
                        pobs->L1[ep-1][sv_G[j]-1]=obs_temp[loc[0]];
                        pobs->L2[ep-1][sv_G[j]-1]=obs_temp[loc[1]];
                        pobs->P1[ep-1][sv_G[j]-1]=obs_temp[loc[2]];
                        pobs->P2[ep-1][sv_G[j]-1]=obs_temp[loc[3]];

                        //free(obs_temp);
                    }
                    // 释放sv_G的内存
                    //free(sv_G);
                }
                //free(loc);
            }
        }
        // 单个文件读取结束，关闭文件
        fclose(file);

        char sav_dat_filename[256]; // 假设文件名长度不超过256个字符

        //生成文件名，路径为当前路径下的RINEX_output_files文件夹
        strcpy(sav_dat_filename, r_opath);
#ifdef LINUX
        strcat(sav_dat_filename, "/");
#elif defined(WINDOWS)
        strcat(sav_dat_filename, "\\");
#else
#endif
        strcat(sav_dat_filename, psitesInfo->name[i]);
        strcat(sav_dat_filename, ".dat");

        writeObsToFile(sav_dat_filename, pobs, 2880, 2880, 2880, 2880, N_SAT);

        //printf("Data %s written to file successfully.\n", sav_dat_filename);

    }
    printf("----------Step 1: completings !----------\n");
}

void read_rinex_v304(const char* r_ipath, const char* r_opath, SitesInfo *psitesInfo, Obs *pobs,const char* OBS_type_1, const char* OBS_type_2, const char* Constellation){    printf("----------Step 1: starting read rinex files!----------\n");
    DIR *dir;
    FILE *file;
    struct dirent *entry;

    // 打开文件夹
    dir = opendir(r_ipath);
    if (dir == NULL) {
        perror("无法打开文件夹");
        exit(EXIT_FAILURE);
    }

    printf("Starting Reading Rinex Files!\n");

    r_file_num = countFilesInDirectory(r_ipath);    // 计算rinex文件数
    char** fileNames;   //存放rinex文件夹下的所有文件名
    malloc_Char_Vector(&fileNames, r_file_num, 256);    //为fileNames分配空间
    getFileNamesInDirectory(r_ipath, fileNames, r_file_num);    //获取rinex文件夹下的所有文件名
    sort_filenames(fileNames, r_file_num);  //对fileNames文件名进行自然排序

    // 为 SitesInfo 结构体内的变量分配内存
    malloc_Char_Vector(&psitesInfo->name, r_file_num, 256);
    malloc_Int_Vector(&psitesInfo->doy, r_file_num);
    malloc_Double_2D(&psitesInfo->coor, r_file_num, 3);
    malloc_Double_Vector(&psitesInfo->RDCB_REF, r_file_num);



    char doy[4];    //用于存放文件名中的doy部分
    char sub_str[3];    //用于存放文件名中的
    // 循环rinex文件名，取出测站名和doy存入psitesInfo里（读sitesInfo）
    for (int i = 0; i < r_file_num; i++) {
        strcpy(psitesInfo->name[i], fileNames[i]);  //fileNames不处理，直接存入name
        //将fileNames[i]的第10到11个字符赋值给sub_str，再将fileNames[i]的第5到7个字符赋值给doy
        //例如，文件名叫gope0010.10o，sub_str就是10，doy就是001
        //现在文件名为：ABMF00GLP_R_20210290000_01D_30S_MO.rnx，sub_str就是21，doy就是029
        strncpy(sub_str, &fileNames[i][14], 2);
        strncpy(doy, &fileNames[i][16], 3);
        sub_str[2] = '\0';  //这一步必须加入字符串结束标志，否则会出问题
        doy[3] = '\0';  //这一步必须加入字符串结束标志，否则会出问题
        psitesInfo->doy[i] = 1000 * atoi(sub_str) + atoi(doy);  //计算真正的doy，存入doy
    }
    // 释放fileNames的内存
    //free(fileNames);
//    //打印psitesInfo中存好的测站名和doy，用于debug
//    for (int i = 0; i < r_file_num; i++) {
//        printf("File No.%d :%s, DOY: %d\n", i + 1, psitesInfo->name[i], psitesInfo->doy[i]);
//    }

    char line[256]; // 假设读文件过程中一行最多包含256个字符
    char filename[256]; //待读的文件名称
    int obst_n; char **obst;
    double h; double m; double s; int ep = 0;
    int *loc;

    // 循环读取每个rinex（读obs观测值）
    for (int i = 0; i < r_file_num; i++) {
        //printf("Reading File No.%d :%s\n", i + 1, psitesInfo->name[i]);

        //拼接部分读取文件的路径
        strcpy(filename, r_ipath);
#ifdef LINUX
        strcat(filename, "/");  //Linux文件系统和windows文件路径符号不一样
#elif defined(WINDOWS)
        strcat(filename, "\\");
#endif
        strcat(filename, psitesInfo->name[i]);  //拼接待读文件完整路径
        file = fopen(filename, "r");    // 打开文件
        if (file == NULL) {
            perror("Can't Open File");
        }

        // 逐行读取
        while (fgets(line, sizeof(line), file)) {
            // -----get receivers coordinates-----------
            if (strstr(line, "APPROX POSITION XYZ") != NULL) {
                double x, y, z;
                // 将行解析成三个双精度浮点数
                if (sscanf(line, "%lf %lf %lf", &x, &y, &z) == 3) {
                    psitesInfo->coor[i][0] = x;
                    psitesInfo->coor[i][1] = y;
                    psitesInfo->coor[i][2] = z;
                }
            }
            //-----get the GPS observables' types------
            // 1.扫描constellation对应的SYS / # / OBS TYPES行，找到对应的OBS_type_1和OBS_type_2对应的index
            // 2.然后开始扫描正文，当是C开头的时候，扫描：1-卫星号；2-OBS_type_1_P1；3-OBS_type_1_L1;4-OBS_type_2_P2;5-OBS_type_2_L2
            // 3.将这些数据存入pobs中
            if (strstr(line, "SYS / # / OBS TYPES") != NULL){
                //扫描第一个字符，如果是constellation，先扫描观测数量，然后建立一个字符串数组存放观测类型
                //扫描每个%3s，一次存入字符串数组
                //扫描下一行，如果第一个字符是空，则继续扫描，直到读取完所有的观测类型
                //然后根据OBS_type_1和OBS_type_2找到对应的index，如果没找到，则报错

                //扫描第一个字符first_char
                char first_char[2];
                strncpy(first_char, &line[0], 1);
                first_char[1] = '\0';

                //如果与Constellation相同，则扫描观测类型数量
                if (strcmp(first_char, Constellation) == 0){
                    //扫描观测数量，第5-6个字符代表观测数量
                    char obst_num[3];
                    strncpy(obst_num, &line[4], 2);
                    obst_num[2] = '\0';
                    obst_n = atoi(obst_num);

                    //建立字符串数组存放观测类型
                    malloc_Char_Vector(&obst, obst_n, 4);
                    int offset = 0;

                    if (obst_n>13){
                        //先读13个
                        offset = 7;
                        for (int j = 0; j < 13; j++){
                            if (sscanf(line + offset, "%3s", obst[j]) == 1) {
                                // 成功解析一个观测名称，更新偏移量，准备解析下一个观测名称
                                offset += 4; // 每个观测名称占据4个字符，包括空格
                            } else {
                                // 解析失败，处理异常情况
                                printf("Failed to parse observation type at index %d\n", j);
                                break;
                            }
                        }
                        //换行，再读剩下的
                        fgets(line, sizeof(line), file);
                        offset = 7;
                        for (int j = 13; j < obst_n; j++){
                            if (sscanf(line + offset, "%3s", obst[j]) == 1) {
                                // 成功解析一个观测名称，更新偏移量，准备解析下一个观测名称
                                offset += 4; // 每个观测名称占据4个字符，包括空格
                            } else {
                                // 解析失败，处理异常情况
                                printf("Failed to parse observation type at index %d\n", j);
                                break;
                            }
                        }
                    }else{  //如果不用换行(<13个），则直接读完就行
                        //开始扫描每个%3s，每次扫描到都存入字符串数组
                        offset = 7;
                        for (int j = 0; j < obst_n; j++){
                            if (sscanf(line + offset, "%3s", obst[j]) == 1) {
                                // 成功解析一个观测名称，更新偏移量，准备解析下一个观测名称
                                offset += 4; // 每个观测名称占据4个字符，包括空格
                            } else {
                                // 解析失败，处理异常情况
                                printf("Failed to parse observation type at index %d\n", j);
                                break;
                            }
                        }
                    }

                    //根据OBS_type_1和OBS_type_2找到对应的index
                    malloc_Int_Vector(&loc, 4);
                    for (int j = 0; j < obst_n; j++){
                        //比较obst[j]和"OBS_type_1"字符串是否相同，如果相同，则给loc[0]赋值为j
                        if (strcmp(obst[j], OBS_type_1) == 0){
                            loc[2] = j;     //P1
                            loc[0] = j+1;   //L1
                        }
                        if (strcmp(obst[j], OBS_type_2) == 0) {
                            loc[3] = j;     //P2
                            loc[1] = j+1;   //L2
                        }
                    }

                    // 释放obst
                    //free(obst);
                }
            }

            //----start get GPS observables------------
            //1.如果第一个字符是">"，则读取时间计算ep
            double buff;
            if (strstr(line, "END OF HEADER") != NULL) {
                while (1) {
                    //读取下一行
                    fgets(line, sizeof(line), file);
                    // 检查是否到达文件结尾
                    if (feof(file)) {
                        // 到达文件结尾，跳出循环
                        break;
                    }
                    //如果第一个字符是">"，则读取时间计算ep
                    if (line[0] == '>') {
                        //将line的第1到2个字符转换成double类型，再赋值给h
                        h = 0;
                        m = 0;
                        s = 0;
                        if (sscanf(line + 13, "%2lf", &h) == 1) {
                        }
                        if (sscanf(line + 16, "%2lf", &m) == 1) {
                        }
                        if (sscanf(line + 19, "%10lf", &s) == 1) {
                        }
                        ep = h * 120 + m * 2 + s / 30 + 1;

                        //---get satellite number----
                        int nsat;
                        int *sv_G;
                        //读取有多少颗卫星观测量，下面就读多少行
                        sscanf(line + 33, "%2d", &nsat);
                        for (int j = 0; j < nsat; j++) {
                            //读取下一行，如果字符与Constellation相同，则开始获取观测值
                            fgets(line, sizeof(line), file);
                            if (line[0] == Constellation[0]) {
                                //读取卫星号
                                int sv;
                                sscanf(line + 1, "%2d", &sv);
                                if (sv==31 || sv>46){
                                    continue;
                                }
                                //读取观测值
                                double* obs_temp;
                                malloc_Double_Vector(&obs_temp, 4);

                                //如果对应位置不是空字符串，则读取对应的观测值
                                for (int k = 0; k<4; k++){
                                    //比较line + 3 + 16 * loc[0], "              ",如果line的第16*loc[0]到第16*loc[0]+14个字符是空格，则给obs_temp[0]赋值为0
                                    if (strncmp(line + 3 + 16 * loc[k], "              ", 14) == 0) {
                                        obs_temp[k] = 0.0;
                                        continue;
                                    } else {
                                        sscanf(line + 3 + 16 * loc[k], "%14lf", &obs_temp[k]);
                                    }
                                }

                                if (sv>31){
                                    pobs->L1[ep - 1][sv - 2] = obs_temp[0];
                                    pobs->L2[ep - 1][sv - 2] = obs_temp[1];
                                    pobs->P1[ep - 1][sv - 2] = obs_temp[2];
                                    pobs->P2[ep - 1][sv - 2] = obs_temp[3];
                                }
                                if (sv<31){
                                    pobs->L1[ep - 1][sv - 1] = obs_temp[0];
                                    pobs->L2[ep - 1][sv - 1] = obs_temp[1];
                                    pobs->P1[ep - 1][sv - 1] = obs_temp[2];
                                    pobs->P2[ep - 1][sv - 1] = obs_temp[3];
                                }

                                free(obs_temp);
                            }
                        }
                    } else {
                        continue;
                    }
                }
            }
        }
        //free(loc);
        // 单个文件读取结束，关闭文件
        fclose(file);

        char sav_dat_filename[256]; // 假设文件名长度不超过256个字符

        //生成文件名，路径为当前路径下的RINEX_output_files文件夹
        strcpy(sav_dat_filename, r_opath);
#ifdef LINUX
        strcat(sav_dat_filename, "/");
#elif defined(WINDOWS)
        strcat(sav_dat_filename, "\\");
#else
#endif
        strcat(sav_dat_filename, psitesInfo->name[i]);
        strcat(sav_dat_filename, ".dat");

        writeObsToFile(sav_dat_filename, pobs, 2880, 2880, 2880, 2880, N_SAT);

        printf("Data %s written to file successfully.\n", sav_dat_filename);

    }
    printf("----------Step 1: completings !----------\n");
}


void read_sp3(const char* s_ipath, const char* s_opath){
    printf("----------Step 2: starting read SP3 files!----------\n");
    DIR *dir;
    FILE *file;
    struct dirent *entry;
    int sp3_file_num;

    // 计算sp3_files文件夹中文件数
    sp3_file_num = countFilesInDirectory(s_ipath);
    //如果文件数小于3，则报错"Need at least 3 SP3 files!"并退出程序
    if (sp3_file_num < 3){
        printf("Need at least 3 SP3 files!\n");
        exit(EXIT_FAILURE);
    }
    char** sp3_filenames;
    malloc_Char_Vector(&sp3_filenames, sp3_file_num, 256);
    getFileNamesInDirectory(s_ipath, sp3_filenames, sp3_file_num);
    sort_filenames(sp3_filenames, sp3_file_num);

    for (int i = 0; i < (sp3_file_num - 2); i++) {
        int index = 0;//读取文件夹中的文件的索引
        char pre_file_name[256];
        char cur_file_name[256];
        char next_file_name[256];
        int G_Week = 0;
        int Day_of_Week = 0;

        //给pre_file_name,cur_file_name,next_file_name都赋值
        for (int j=0; j<sp3_file_num; j++) {  // 遍历文件夹中的每一个文件
            if (index == i) {  //前两个文件分别是"."和"..",所以从第三个文件开始读
                char filepath[256];  // 定义一个足够大的缓冲区来存储完整路径
                strcpy(filepath, s_ipath);  // 将文件夹路径拷贝到缓冲区

#ifdef LINUX
                strcat(filepath, "/");
#elif defined(WINDOWS)
                strcat(filepath, "\\");  // 将文件名拷贝到缓冲区
#else
#endif
                strcat(filepath, sp3_filenames[j]);  // 将文件名拷贝到缓冲区
                //把filepath的值赋给pre_file_name
                strcpy(pre_file_name, filepath);
            }
            if (index == i+1){
                char filepath[256];  // 定义一个足够大的缓冲区来存储完整路径
                strcpy(filepath, s_ipath);  // 将文件夹路径拷贝到缓冲区

#ifdef LINUX
                strcat(filepath, "/");
#elif defined(WINDOWS)
                strcat(filepath, "\\");  // 将文件名拷贝到缓冲区
#else
#endif
                strcat(filepath, sp3_filenames[j]);  // 将文件名拷贝到缓冲区
                strcpy(cur_file_name, filepath);

                sscanf(sp3_filenames[j] + 3, "%4d", &G_Week);
                sscanf(sp3_filenames[j] + 7, "%1d", &Day_of_Week);
            }
            if (index == i+2){
                char filepath[256];  // 定义一个足够大的缓冲区来存储完整路径
                strcpy(filepath, s_ipath);  // 将文件夹路径拷贝到缓冲区

#ifdef LINUX
                strcat(filepath, "/");
#elif defined(WINDOWS)
                strcat(filepath, "\\");  // 将文件名拷贝到缓冲区
#else
#endif
                strcat(filepath, sp3_filenames[j]);  // 将文件名拷贝到缓冲区
                strcpy(next_file_name, filepath);
            }
                index++;  // 如果不是第 i 个文件，继续查找下一个文件
        }

        //读取pre_file_name sp3文件的x,y,z坐标
        //定义一个一维double数组pre_xyz，长度为3，并赋初值为0
        double ***pre_xyz; double ***cur_xyz; double ***next_xyz; double ***sate_xyz;
        malloc_double_3D(&pre_xyz, 3, 96, N_SAT);
        malloc_double_3D(&cur_xyz, 3, 96, N_SAT);
        malloc_double_3D(&next_xyz, 3, 96, N_SAT);
        malloc_double_3D(&sate_xyz, 3, 2880, N_SAT);

        //pre_xyz[][][],第一个[]代表x,y,z坐标，第二个[]代表行（时间），第三个[]代表列（卫星编号）
        parse_sp3(pre_file_name, pre_xyz);
        parse_sp3(cur_file_name, cur_xyz);
        parse_sp3(next_file_name, next_xyz);

        interplotation(pre_xyz, cur_xyz, next_xyz, sate_xyz);

        int* DOY = GWeek_2_DOY(G_Week, Day_of_Week);

        char sav_filename[256]; // 假设文件名长度不超过256个字符
        //生成文件名，路径为当前路径下的SP3_output_files文件夹
        strcpy(sav_filename, s_opath);
#ifdef LINUX
        strcat(sav_filename, "/");
#elif defined(WINDOWS)
        strcat(sav_filename, "\\");
#else
#endif
        //添加G_Week和Day_of_Week到文件名中
        char G_Week_str[5];
        char Day_of_Week_str[2];
        sprintf(G_Week_str, "%d", DOY[0]);
        sprintf(Day_of_Week_str, "%d", DOY[1]);
        strcat(sav_filename, G_Week_str);
        strcat(sav_filename, "_");
        strcat(sav_filename, Day_of_Week_str);
        strcat(sav_filename, "sp3");
        strcat(sav_filename, ".dat");

        //saveSp3ToCSV(sate_xyz, sav_filename); //以csv格式存储sp3文件
        save_sp3_To_Bin_File(sate_xyz, 3, 2880, N_SAT, sav_filename);  //以二进制格式存储sp3文件

        free(pre_xyz);
        free(cur_xyz);
        free(next_xyz);
        free(sate_xyz);
        free(DOY);
    }
    free(sp3_filenames);
    printf("----------Step 2: completings !----------\n");
}
// free到这儿了
// 读取ionex文件，相关结果存放在sitesInfo.RDCB_REF和sDCB_REF结构体中
void read_ionex(const char* i_ipath, const char* i_opath, SitesInfo *psitesInfo, SDCB_REF sdcb_ref){
    printf("----------Step 3: starting read ionex files!----------\n");
    int n = r_file_num;
    int new_size = 0;
    int* doys; //去重的DOY存放在doys数组中，大小为new_size
    //给doys申请n个int类型的内存，并赋初值为0
    malloc_Int_Vector(&doys, n);

    //遍历r_file的doy
    doys[0] = psitesInfo->doy[0];
    new_size++;
    for (int i = 1; i < n; i++){
        //查看doys数组中是否有与psitesInfo->doy[i]相同的元素，若有，则不做任何行为，若无，则将psitesInfo->doy[i]赋值给doys[new_size]，并且new_size加1
        int flag = 0;
        for (int j = 0; j < new_size; j++){
            if (doys[j] == psitesInfo->doy[i]){
                flag = 1;   //说明在这个小循环中找到了相同的元素
                break;
            }
        }
        if (flag == 0){
            doys[new_size] = psitesInfo->doy[i];
            new_size++;
        }
    }

    //为sdcb_ref结构体内的数组分配内存,value需要((new_size)*32)的尺寸
    malloc_Double_2D(&sdcb_ref.value, new_size, N_SAT);
    malloc_Int_Vector(&sdcb_ref.doy, new_size);

    char** ionex_filenames; int ionex_file_num = 0;

    ionex_file_num = countFilesInDirectory(i_ipath);

    //给ionex_filenames申请ionex_file_num个char*类型的内存,每个char不超过32个字符
    malloc_Char_Vector(&ionex_filenames, ionex_file_num, N_SAT);

    getFileNamesInDirectory(i_ipath, ionex_filenames, ionex_file_num);
    sort_filenames(ionex_filenames, ionex_file_num);

    //为DCB_rec分配内存, 大小为ionex_file_num
    malloc_Double_Vector(&DCB_rec, ionex_file_num);

    //开始读文件
    for (int i= 0; i < new_size; i++){
        int index = 0; int *index2; int doy_num;
        char** sites_name;

        malloc_Int_Vector(&index2, r_file_num);

        index = find_ionex_index(doys[i], ionex_filenames, ionex_file_num);
        doy_num = find_rinex_indexs(index2, doys[i]); //在rinex文件列表中找到所有和doy[i]相等的索引号, doy_num为索引到相同值的个数

        //给sites_name申请doy_num个char*类型的内存,每个char不超过4个字符
        malloc_Char_Vector(&sites_name, doy_num, 4);

        //取出sites_name的值
        for (int j = 0; j < doy_num; j++){
            strncpy(sites_name[j], &psitesInfo->name[index2[j]][0],4);
            sites_name[j][4] = '\0';
        }

        //开始读ionex文件
        char ionex_file_path[256];
        strcpy(ionex_file_path, i_ipath);
#ifdef LINUX
        strcat(ionex_file_path, "/");
#elif defined(WINDOWS)
        strcat(ionex_file_path, "\\");
#else
#endif
        strcat(ionex_file_path, ionex_filenames[index]);

        read_single_ionex(ionex_file_path, i, sites_name, doy_num, sdcb_ref, DCB_rec);

        for (int j = 0; j < ionex_file_num; j++){
            sitesInfo.RDCB_REF[index2[j]] = DCB_rec[j];
        }

        sdcb_ref.doy[i] = doys[i];

        free(index2);
        free(sites_name);
    }
    free(DCB_rec);
    free(ionex_filenames);
    free(doys);
    printf("----------Step 3: completings !----------\n");
}

// 文件夹内文件数
int countFilesInDirectory(const char *folderPath) {
    DIR *dir;
    struct dirent *entry;
    int counter = 0;

    // 打开文件夹
    dir = opendir(folderPath);

    if (dir == NULL) {
        perror("无法打开文件夹");
        exit(EXIT_FAILURE);
    }

    // 读取文件夹中的文件
    while ((entry = readdir(dir)) != NULL) {
        if (entry->d_name[0] != '.'){
            counter ++;
            //printf("No.%d :%s\n", counter,entry->d_name); //Print file names
        }
    }

    // 关闭文件夹
    closedir(dir);

    return counter;
}

// 功能描述：只适用于将obs文件写入csv文件
// 入参：*file —— 存到什么文件路径下的什么文件； *array —— 待存储的二维数组；
//      * dataType —— 数据标签，自己定义；
// 出参：无
void writeArrayToFileWithLabels(FILE *file, double **array, const char* dataType) {
    // 写入列标签
    fprintf(file, "%s_Epoch", dataType); // 假设第一列为时间戳
    for (int j = 1; j <= N_SAT; j++) {
        fprintf(file, ",%s_sv_%d", dataType, j); // 生成并写入标签
    }
    fprintf(file, "\n"); // 结束标签行

    // 写入数据
    for (int i = 0; i < 2880; i++) {
        fprintf(file, "%d", i+1); // 假设第一列为时间戳，此处简化处理
        for (int j = 0; j < N_SAT; j++) {
            fprintf(file, ",%.4f", array[i][j]);
        }
        fprintf(file, "\n");
    }
}

// 功能描述：解析sp3文件，将x、y、z上的2880*32都存入xyz三维数组，每个维度分别代表x、y、z
// 入参：* sp3_file —— sp3文件完整路径； ***xyz —— 待存放精密星历xyz坐标的三维数组
// 出参：无
void parse_sp3(char* sp3_file, double ***xyz){
    FILE *file;
    char line[256];
    file = fopen(sp3_file, "r");
    if (file == NULL) {
        perror("Can't Open File");
    }

    int ep = 0;

    while (fgets(line, sizeof(line), file)) {
        //如果line开头第一个字符是"*"，则执行下面的语句
        if (line[0] == '*'){
            //读取line的第15到16个字符，转换为double类型，赋给一个变量h
            double h = 0; double m = 0;
            if (sscanf(line + 14, "%2lf", &h) == 1){
            }
            if (sscanf(line + 17, "%2lf", &m) == 1){
            }
            //对m/15的结果进行四舍五入
            ep = (int)h*4 + round(m/15) + 1;
            continue;
        }
        //如果line长度大于1，并且line的第一个和第二个字符是"PG"，则执行下面的语句
        if (strlen(line) > 1 && line[0] == 'P' && line[1] == 'G'){
            int sv = 0;
            if (sscanf(line+2, "%2d", &sv) == 1){
            }
            if (sscanf(line+4, "%14lf", &xyz[0][ep-1][sv-1]) == 1){
            }
            if (sscanf(line+18, "%14lf", &xyz[1][ep-1][sv-1]) == 1){
            }
            if (sscanf(line+32, "%14lf", &xyz[2][ep-1][sv-1]) == 1){
            }
            continue;
        }
    }
    fclose(file);
}

// 功能描述：插值，一天的精密坐标每15min一组数据，一天的数据量就是(3*96*32)，结合前一天、后一天的精密坐标，
//         插值生成(3*2880*32)尺寸的精密坐标（每30s一组数据）。目的是与rinex观测值数据量对齐
// 入参：***pre_xyz —— 上一时刻的xyz坐标； ***cur_xyz —— 当前时刻的xyz坐标； ***next_xyz —— 下一时刻的xyz坐标；
//      ***sate_xyz —— 插值后的xyz坐标
// 出参：无
void interplotation(double ***pre_xyz, double ***cur_xyz, double ***next_xyz, double ***sate_xyz){
    double ***cur_xyz_etd;  //当前xyz的扩展数组
    malloc_double_3D(&cur_xyz_etd, 3, 105, N_SAT); //cur_xyz_etd的尺寸为3*105*32

    extend_matrix(0, pre_xyz, cur_xyz, next_xyz, cur_xyz_etd);//x_extend
    extend_matrix(1, pre_xyz, cur_xyz, next_xyz, cur_xyz_etd);//y_extend
    extend_matrix(2, pre_xyz, cur_xyz, next_xyz, cur_xyz_etd);//z_extend

    //创建一个整数数组m_t，起始值为-120，结束值为3000，数组长度为105
    int m_t[105];
    for (int i = 0; i < 105; i++){
        m_t[i] = -120 + 30*i;
    }

    for (int i = 0; i < N_SAT; i++){
        for (int j = 0; j < 96; j++){
            double tt[10];
            //tt的值等于m_t的第j个元素到第j+9个元素
            for (int k = 0; k < 10; k++){
                tt[k] = m_t[j+k];
            }
            //创建一个int类型的数组x_temp，cur_xyz_etd[0]的第j行到第j+9行，第i列的值赋给x_temp
            double x_temp[10];
            for (int k = 0; k < 10; k++){
                x_temp[k] = cur_xyz_etd[0][j+k][i];
            }
            //创建一个int类型的数组y_temp，cur_xyz_etd[1]的第j行到第j+9行，第i列的值赋给y_temp
            double y_temp[10];
            for (int k = 0; k < 10; k++){
                y_temp[k] = cur_xyz_etd[1][j+k][i];
            }
            //创建一个int类型的数组z_temp，cur_xyz_etd[2]的第j行到第j+9行，第i列的值赋给z_temp
            double z_temp[10];
            for (int k = 0; k < 10; k++){
                z_temp[k] = cur_xyz_etd[2][j+k][i];
            }
            //创建一个t0数组，起始值为m_t(j+4),终点值为m_t(j+5)-1，长度为30
            double t0[30];
            for (int k = 0; k < 30; k++){
                t0[k] = m_t[j+4] + k;
            }

            double* x0 = interp_lag(tt, x_temp, t0);
            double* y0 = interp_lag(tt, y_temp, t0);
            double* z0 = interp_lag(tt, z_temp, t0);

            for (int k = 0; k < 30; k++){
                sate_xyz[0][30*(j+1)-30+k][i] = x0[k];
                sate_xyz[1][30*(j+1)-30+k][i] = y0[k];
                sate_xyz[2][30*(j+1)-30+k][i] = z0[k];
            }

            free(x0);
            free(y0);
            free(z0);
        }
    }

    free(cur_xyz_etd);

}

// 功能描述：获得一个扩展矩阵，从矩阵pre_xyz[0]中选取第93到第96行的所有列，形成一个子矩阵。将子矩阵与矩阵cur_xyz[0]进行垂直拼接，
//         即将子矩阵放在cur_xyz[0]的上方。将矩阵next_xyz[0]中的第1到第5行的所有列与前一步得到的结果进行垂直拼接，即将next_xyz[0]
//         的前5行放在前一步得到的结果的下方。
// 入参：***pre_xyz —— 上一时刻的xyz坐标； ***cur_xyz —— 当前时刻的xyz坐标； ***next_xyz —— 下一时刻的xyz坐标；
//      ***xyz_etd —— 扩展后的矩阵
// 出参：无
void extend_matrix(int dimension, double ***pre_xyz, double ***cur_xyz, double ***next_xyz, double ***xyz_etd){
    //先填充pre_xyz的部分
    for (int j = 0; j < 4; j++){
        for (int k = 0; k < N_SAT; k++) {
            xyz_etd[dimension][j][k] = pre_xyz[dimension][j+92][k];
        }
    }

    // 填充cur_xyz的部分
    for (int j = 4; j < 100; j++) {
        for (int k = 0; k < N_SAT; k++) {
            xyz_etd[dimension][j][k] = cur_xyz[dimension][j-4][k];
        }
    }

    // 填充next_xyz的部分
    for (int j = 100; j < 105; j++) {
        for (int k = 0; k < N_SAT; k++) {
            xyz_etd[dimension][j][k] = next_xyz[dimension][j-100][k];
        }
    }
}

// 功能描述：
// 入参：
// 出参：无
double* interp_lag(double* x, double* y, double* x0) {
    int n = 10; // 假设x中有10个点
    int n0 = 30;
    double* y0;
    malloc_Double_Vector(&y0, n0);
    if (y0 == NULL) { // 检查内存分配是否成功
        printf("Memory allocation failed.\n");
        return NULL;
    }

    for (int k = 0; k < n0; ++k) { // 遍历x0中的每个点
        y0[k] = 0; // 初始化y0的当前元素
        for (int i = 0; i < n; ++i) { // 遍历x中的每个点
            double t = 1;
            for (int j = 0; j < n; ++j) { // 计算插值基函数
                if (j != i) {
                    t *= (x0[k] - x[j]) / (x[i] - x[j]);
                }
            }
            y0[k] += t * y[i]; // 累加计算插值结果
        }
    }

    return y0; // 返回计算结果
}

// 功能描述：通过卫星周和DOW计算当前DOY
// 入参：G_Week —— 卫星周（GPS周）； Day_of_Week —— 一个卫星周中的第X天
// 出参：int DOY —— 一年中的第x天
int* GWeek_2_DOY(int G_Week, int Day_of_Week){
    int *DOY;
    malloc_Int_Vector(&DOY, 2);
    int doy = 0;
    int year = 1980;
    int n=0;

    int a[] = {360,725,1090,1455,1821,2186,2551,2916,3282,3647,4012,4377,4743,5108,5473,5838,6204,6569,6934,7299,7665,8030,8395,8760,9126,9491,9856,10221,10587,10952,11317,11682,12048,12413,12778,13143,13509,13874,14239,14604,14970};
    n = G_Week*7 + Day_of_Week;

    for (int i = 0; i < 41; i++){
        if (!(n > a[0])){
            doy = n + 6;
            break;
        }
        if (!(n > a[i])){
            year = year + i;
            doy = n - a[i-1];
            break;
        }
    }

    DOY[0] = year;
    DOY[1] = doy;
    return DOY;
}

void saveSp3ToCSV(double*** sate_xyz, char* filename) {
    FILE* file = fopen(filename, "w"); // 打开文件用于写入
    if (file == NULL) {
        printf("无法打开文件 %s\n", filename);
        return;
    }

    // 为卫星号SVN创建表头
    fprintf(file, "Epoch/SVN,");
    for (int svn = 1; svn <= N_SAT; svn++) {
        fprintf(file, "%d", svn);
        if (svn < N_SAT) fprintf(file, ",");
        else fprintf(file, "\n");
    }

    // 分别按照x、y、z的顺序写入数据
    for (int dim = 0; dim < 3; dim++) { // 0: x, 1: y, 2: z
        for (int epoch = 0; epoch < 2880; epoch++) {
            fprintf(file, "Epoch %d,", epoch + 1);
            for (int svn = 0; svn < N_SAT; svn++) {
                fprintf(file, "%lf", sate_xyz[dim][epoch][svn]);
                if (svn < 31) fprintf(file, ",");
                else fprintf(file, "\n");
            }
        }
        // 在不同维度数据间添加空行，增加可读性，可根据个人喜好决定是否保留
        fprintf(file, "\n");
    }

    fclose(file); // 关闭文件
    //printf("SP3 data successfully saved to %s\n", filename);
}

int find_ionex_index(int doy, char** ionex_filenames, int ionex_file_num){
    int index = 0;

    //读取i_ipath文件夹中的文件
    for(int i = 0; i < ionex_file_num; i++){
        int number1 = 0; int number2 = 0; int doy_temp = 0;

        //第5到第7位的int类型数字1
        sscanf(ionex_filenames[i] + 4, "%3d", &number1);
        //第10到第11位的int类型数字2
        sscanf(ionex_filenames[i] + 9, "%2d", &number2);
        doy_temp = number1 + number2*1000;
        if (doy_temp == doy){
            return index;
            break;
        }
        index++;
    }
    return index;
}

int find_rinex_indexs(int* index, int doy){
    int counter = 0;

    for(int i = 0; i < r_file_num; i++){
        if (sitesInfo.doy[i] == doy){
            index[counter] = i;
            counter++;
        }
    }

    return counter;
}

//读取单个ionex文件，将sDCB和rDCB存入sDCB_REF.value和DCB_rec中
void read_single_ionex(const char* ionex_file, int i, char** sites_name, int doy_num, SDCB_REF sDCB_REF, double* DCB_rec) {
    int flag = 0;
    int prn = 0;

    // 打开文件
    FILE *file;
    file = fopen(ionex_file, "r");
    if (file == NULL) {
        perror("Can't Open File");
    }

    // 读取文件
    char line[256];
    while (fgets(line, sizeof(line), file)) {
        if (flag == 1){
            break;
        }

        //如果line长度大于76并且第61到77个字符是"START OF AUX DATA"，则执行下面的语句
        if (strlen(line) > 76 && strncmp(line + 60, "START OF AUX DATA", 17) == 0){
            flag = 1;
            while (fgets(line, sizeof(line), file)) {
                //如果line长度大于74并且第61到75个字符是"END OF AUX DATA"，则跳出循环
                if (strlen(line) > 74 && strncmp(line + 60, "END OF AUX DATA", 15) == 0){
                    break;
                }
                //--satellites' DCB
                //如果line长度大于75，并且第4个字符是"G"或空格，并且第61到75个字符是"PRN / BIAS / RMS"，则执行下面的语句
                if (strlen(line) > 75 && (line[3] == 'G' || line[3] == ' ') && !(strncmp(line + 60, "PRN / BIAS / RMS", 15)) == 1){
                    sscanf(line + 4, "%2d", &prn);
                    sscanf(line + 9, "%7lf", &sDCB_REF.value[i][prn-1]);
                    continue;
                }
                //--receivers' DCB
                if (strlen(line) > 10 && (line[3] == 'G' || line[3] == ' ') && !(strncmp(line + 60, "STATION / BIAS / RMS", 20)) == 1){
                    //比较line的第7到第10个字符是否在sites_name中，如果在，则执行下面的语句
                    for (int j = 0; j < doy_num; j++){
                        if (strncmp(line + 6, sites_name[j], 4) == 0){
                            sscanf(line + 29, "%7lf", &DCB_rec[j]);
                            break;
                        }
                    }
                }
            }
        }
    }
    //关闭文件
    fclose(file);
}

void get_smoothed_P4(SitesInfo sitesInfo, double z_threshold, int flag){
    printf("----------Step 4: starting get smoothed P4!----------\n");
    //找到去重后的rinex文件列表中与site对应的index
    int new_size = 0;
    char** sites; //去重的DOY存放在doys数组中，大小为new_size

    //通过一个函数得到一个只包含sitesInfo.name中前四位（site名称）的数组siteOnlyName，大小和sitesInfo.name一样
    char** sitesOnlyName;
    malloc_Char_Vector(&sitesOnlyName, r_file_num, 4);

    for (int i = 0; i < r_file_num; i++){
        strncpy(sitesOnlyName[i], &sitesInfo.name[i][0],4);
        sitesOnlyName[i][4] = '\0';
    }

    removeSameElementAndReturn_Char(sitesOnlyName, &sites, r_file_num, &new_size);
    free(sitesOnlyName);

    for (int j = 1; j < r_file_num; j++){
        //查看doys数组中是否有与psitesInfo->doy[i]相同的元素，若有，则不做任何行为，若无，则将psitesInfo->doy[i]赋值给doys[new_size]，并且new_size加1
        int flag = 0;
        char* temp_name;
        temp_name = (char *) malloc(4 * sizeof(char));
        sscanf(sitesInfo.name[j], "%4s", temp_name);

        for (int k = 0; k < new_size; k++){
            //如果sites[j]和temp_name字符串完全相同
            if (strcmp(sites[k], temp_name) == 0){
                flag = 1;   //说明在这个小循环中找到了相同的元素
                break;
            }
        }
        if (flag == 0){
            sscanf(temp_name, "%4s", sites[new_size]);
            new_size++;
        }

        free(temp_name);
    }

    //依次读取每天的rinex参数和sp3参数
    for (int i=0; i < r_file_num; i++){
        //读取RINEX_output_files文件夹中的第一个文件
        char rinex_file_path[256];
        Obs obs_temp;
        strcpy(rinex_file_path, r_opath);
#ifdef LINUX
        strcat(rinex_file_path, "/");
#elif defined(WINDOWS)
        strcat(rinex_file_path, "\\");
#else
#endif
        strcat(rinex_file_path, sitesInfo.name[i]);
        strcat(rinex_file_path, ".dat");
        readObsFromFile(rinex_file_path, &obs_temp);

        char* site; int doy;
        site = (char *) malloc(64 * sizeof(char));
        sscanf(sitesInfo.name[i], "%4s", site);
        doy = sitesInfo.doy[i];
        int index = 0;

        for (int j=0;j<new_size;j++){
            if (strcmp(sites[j], site) == 0){
                index = j;
                break;
            }
        }

        double sx = 0; double sy = 0; double sz = 0;
        sx = sitesInfo.coor[index][0];
        sy = sitesInfo.coor[index][1];
        sz = sitesInfo.coor[index][2];

        if (sx == 0 && sy == 0 && sz == 0){
            continue;
        }

        //组成读取sp3的文件名
        char sp3_load_file_path[256];

        //添加s_opath到sp3_load_file_path
        strcpy(sp3_load_file_path, s_opath);
#ifdef LINUX
        strcat(sp3_load_file_path, "/");
#elif defined(WINDOWS)
        strcat(sp3_load_file_path, "\\");
#else
#endif
        //在字符串中添加"20"
        strcat(sp3_load_file_path, "20");
        //在字符串sp3_load_file_path中添加doy的前2个字符，如doy=10001，则添加“10”
        // 将doy的前两个字符转换为字符串
        char doy_str[3];
        snprintf(doy_str, sizeof(doy_str), "%d", doy);
        strncat(sp3_load_file_path, doy_str, 2);
        //在字符串中添加"_"
        strcat(sp3_load_file_path, "_");
        // 提取doy的后三位数字
        int last_three_digits = doy % 1000;
        // 转换后三位数字为字符串
        char doy_str_2[4]; // 考虑到最大值为999，因此最多需要4个字符（包括终止符'\0'）
        sprintf(doy_str_2, "%d", last_three_digits);
        strcat(sp3_load_file_path, doy_str_2);
        //在字符串中添加"sp3.dat"
        strcat(sp3_load_file_path, "sp3.dat");

        //读取sp3文件
        double ***sate_xyz;
        malloc_double_3D(&sate_xyz, 3, 2880, N_SAT);

        load_Sp3_From_Bin_File(&sate_xyz, 3, 2880, N_SAT, sp3_load_file_path);

        cutobs(&obs_temp, sate_xyz, sx, sy, sz, z_threshold);//筛选后的obs值存在obs_temp中

        double** P4;
        malloc_Double_2D(&P4, 2880, N_SAT);

        P4 = pre_pro(&obs_temp);//对obs_temp进行预处理

        char sav_P4_filename[256];
        //添加m_p4_path到sav_P4_filename
        strcpy(sav_P4_filename, m_p4_path);
#ifdef LINUX
        strcat(sav_P4_filename, "/");
#elif defined(WINDOWS)
        //添加“\\”到sav_P4_filename
        strcat(sav_P4_filename, "\\");
#else
#endif
        //先将site的值赋给sav_P4_filename
        strcat(sav_P4_filename, site);
        //再取doy赋给sav_P4_filename
        char doy_str_3[6];
        snprintf(doy_str_3, sizeof(doy_str_3), "%d", doy);
        strncat(sav_P4_filename, doy_str_3, 5);
        //在字符串中添加"_P4.dat"
        strcat(sav_P4_filename, "_P4.dat");

        save_P4_To_Bin_File(P4, 2880, N_SAT, sav_P4_filename);//将P4存入文件中

        free(site);
        free(sate_xyz);
        free(P4);
    }
    printf("----------Step 4: completings !----------\n");
}

void writeObsToFile(const char* filename, Obs* obs, int p1Rows, int p2Rows, int l1Rows, int l2Rows, int cols) {
    FILE* file = fopen(filename, "wb");
    if (file == NULL) {
        perror("Failed to open file");
        return;
    }

    // 写入每个数组的行数和列数
    fwrite(&p1Rows, sizeof(int), 1, file);
    fwrite(&p2Rows, sizeof(int), 1, file);
    fwrite(&l1Rows, sizeof(int), 1, file);
    fwrite(&l2Rows, sizeof(int), 1, file);
    fwrite(&cols, sizeof(int), 1, file);

    // 为了简化，这里假设所有数组的列数相同
    for (int i = 0; i < p1Rows; i++) {
        fwrite(obs->P1[i], sizeof(double), cols, file);
    }
    // 重复上述过程写入P2, L1, L2
    for (int i = 0; i < p2Rows; i++) {
        fwrite(obs->P2[i], sizeof(double), cols, file);
    }
    for (int i = 0; i < l1Rows; i++) {
        fwrite(obs->L1[i], sizeof(double), cols, file);
    }
    for (int i = 0; i < l2Rows; i++) {
        fwrite(obs->L2[i], sizeof(double), cols, file);
    }

    fclose(file);
}

void readObsFromFile(const char* filename, Obs* obs) {
    FILE* file = fopen(filename, "rb");
    if (file == NULL) {
        perror("Failed to open file");
        return;
    }

    int p1Rows, p2Rows, l1Rows, l2Rows, cols;
    fread(&p1Rows, sizeof(int), 1, file);
    fread(&p2Rows, sizeof(int), 1, file);
    fread(&l1Rows, sizeof(int), 1, file);
    fread(&l2Rows, sizeof(int), 1, file);
    fread(&cols, sizeof(int), 1, file);

    // 根据读取的行列数分配内存和读取数据
    // 注意处理内存分配失败的情况
    // 以下代码为伪代码，需要根据实际情况调整
    malloc_Double_2D(&(obs->P1), p1Rows, cols);
    malloc_Double_2D(&(obs->P2), p2Rows, cols);
    malloc_Double_2D(&(obs->L1), l1Rows, cols);
    malloc_Double_2D(&(obs->L2), l2Rows, cols);

    // 为每个数组分配内存并从文件中读取数据
    for (int i = 0; i < p1Rows; i++) {
        fread(obs->P1[i], sizeof(double), cols, file);
    }
    for (int i = 0; i < p2Rows; i++) {
        fread(obs->P2[i], sizeof(double), cols, file);
    }
    for (int i = 0; i < l1Rows; i++) {
        fread(obs->L1[i], sizeof(double), cols, file);
    }
    for (int i = 0; i < l2Rows; i++) {
        fread(obs->L2[i], sizeof(double), cols, file);
    }

    fclose(file);
}

void save_sp3_To_Bin_File(double ***sate_xyz, int dim1, int dim2, int dim3, const char *filename) {
    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
        perror("Failed to open file");
        exit(1);
    }

    for (int i = 0; i < dim1; ++i) {
        for (int j = 0; j < dim2; ++j) {
            fwrite(sate_xyz[i][j], sizeof(double), dim3, file);
        }
    }

    fclose(file);
}

void load_Sp3_From_Bin_File(double ****sate_xyz, int dim1, int dim2, int dim3, const char *filename) {
    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
        perror("Failed to open file");
        exit(1);
    }

    *sate_xyz = (double ***)malloc(dim1 * sizeof(double **));
    for (int i = 0; i < dim1; ++i) {
        (*sate_xyz)[i] = (double **)malloc(dim2 * sizeof(double *));
        for (int j = 0; j < dim2; ++j) {
            (*sate_xyz)[i][j] = (double *)malloc(dim3 * sizeof(double));
            fread((*sate_xyz)[i][j], sizeof(double), dim3, file);
        }
    }

    fclose(file);
}

void cutobs(Obs *obs_temp, double*** sate_xyz, double sx, double sy, double sz, double lim){
    double** x = sate_xyz[0];
    double** y = sate_xyz[1];
    double** z = sate_xyz[2];

    for (int i = 0; i < N_SAT; i++){
        for (int j = 0; j < 2880; j++){
            if (obs_temp->L1[j][i] == 0 || obs_temp->L2[j][i] == 0 || obs_temp->P1[j][i] == 0 || obs_temp->P2[j][i] == 0){
                obs_temp->L1[j][i] = 0; obs_temp->L2[j][i] = 0; obs_temp->P1[j][i] = 0; obs_temp->P2[j][i] = 0;
                continue;
            }

            //计算E、A
            double* EA = get_EA(sx,sy,sz,x[j][i]*1000,y[j][i]*1000,z[j][i]*1000);
            double el = EA[0];
            double aaa = EA[1];

            if (el<lim){
                obs_temp->L1[j][i] = 0; obs_temp->L2[j][i] = 0; obs_temp->P1[j][i] = 0; obs_temp->P2[j][i] = 0;
                continue;
            }
        }
    }
}

double* get_EA(double sx, double sy, double sz, double x, double y, double z){
    double* EA; malloc_Double_Vector(&EA, 2);
    double el = 0; double aaa = 0; double SB = 0; double SL = 0;
    double deta_xyz[3]; double* BLH;
    double NEU[3] = {0, 0, 0};

    BLH = XYZ2BLH(sx, sy, sz);
    SB = BLH[0];
    SL = BLH[1];

    //转换矩阵
    //    T=[-sin(sb)*cos(sl) -sin(sb)*sin(sl) cos(sb);
    //    -sin(sl)               cos(sl)         0;
    //    cos(sb)*cos(sl) cos(sb)*sin(sl)  sin(sb)];%transition matrix(XYZ to NEU)
    double T[3][3] = {{-sin(SB)*cos(SL), -sin(SB)*sin(SL), cos(SB)},
                      {-sin(SL), cos(SL), 0},
                      {cos(SB)*cos(SL), cos(SB)*sin(SL), sin(SB)}};
    //deta_xyz=[x,y,z]-[sx,sy,sz];
    deta_xyz[0] = x - sx;
    deta_xyz[1] = y - sy;
    deta_xyz[2] = z - sz;
    //NEU=T*(deta_xyz)';
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            NEU[i] += T[i][j] * deta_xyz[j];
        }
    }
    //E=atan(NEU(3)/sqrt(NEU(1)*NEU(1)+NEU(2)*NEU(2)));
    el = atan(NEU[2]/sqrt(NEU[0]*NEU[0]+NEU[1]*NEU[1]));
    //A=atan(abs(NEU(2)/NEU(1)));
    aaa = atan(fabs(NEU[1]/NEU[0]));
    //
    if (NEU[0]>0){
        if (NEU[1]>0){}
        else{
            aaa = 2*pi - aaa;
        }
    }else{
        if (NEU[1]>0){
            aaa = pi - aaa;
        }else{
            aaa = pi + aaa;
        }
    }

    EA[0] = el;
    EA[1] = aaa;
    return EA;
}

double* XYZ2BLH(double x, double y, double z){
    double* BLH; malloc_Double_Vector(&BLH, 3);
    double B = 0; double B0 = 0; double L = 0; double H = 0; double N = 0;

    double a = 6378137.0;
    double e2=0.0066943799013;
    L=atan(fabs(y/x));

    if (y>0){
        if (x>0){
        }else{
            L= pi - L;
        }
    }else{
        if (x>0){
            L=2*pi-L;
        }else{
            L=pi+L;
        }
    }

    B0 = atan(z/sqrt(x*x+y*y));

    while(1){
        N = a/sqrt(1-e2*sin(B0)*sin(B0));
        H = z/sin(B0)-N*(1-e2);
        B = atan(z*(N+H)/(sqrt(x*x+y*y)*(N*(1-e2)+H)));
        if (fabs(B-B0)<1e-6) {
            break;
        }
        B0 = B;
    }

    N = a/sqrt(1-e2*sin(B)*sin(B));

    BLH[0] = B;
    BLH[1] = L;
    BLH[2] = H;
    return BLH;
}

double** pre_pro(Obs* obs_temp){

    double** L6; double** Li; double** Nw; double** P4; double** L4;
    malloc_Double_2D(&L6, 2880, N_SAT);
    malloc_Double_2D(&Li, 2880, N_SAT);
    malloc_Double_2D(&Nw, 2880, N_SAT);
    malloc_Double_2D(&P4, 2880, N_SAT);
    malloc_Double_2D(&L4, 2880, N_SAT);

    // Wide Lane Wavelength
    double lambda_w = c/(f1-f2);
    // MW observable
    //L6 = lambda_w*(obs_temp->L1 - obs_temp->L2) - (f1*obs.P1+f2*obs.P2)/(f1+f2);
    for (int i = 0; i < 2880; i++){
        for (int j = 0; j < N_SAT; j++){
            L6[i][j] = lambda_w*(obs_temp->L1[i][j] - obs_temp->L2[i][j]) - (f1*obs_temp->P1[i][j]+f2*obs_temp->P2[i][j])/(f1+f2);
        }
    }
    //Li = obs.L1 - f1*obs.L2/f2;
    for (int i = 0; i < 2880; i++){
        for (int j = 0; j < N_SAT; j++){
            Li[i][j] = obs_temp->L1[i][j] - f1*obs_temp->L2[i][j]/f2;
        }
    }

    //Nw = L6/lambda_w;
    for (int i = 0; i < 2880; i++){
        for (int j = 0; j < N_SAT; j++){
            Nw[i][j] = L6[i][j]/lambda_w;
        }
    }

    for (int prn = 0; prn < N_SAT; prn++){
        int** arc; int arc_n = 0 ; int aaa = 2;

        //------divide arc---------------------------
        arc = Get_arc(L6, prn, &arc_n);//arc是一个n行2列的数组，每一行代表一个arc的起始和终止的index

        //----delete arc less than 10 epoches-------
        int* arc_d; int num_of_arc_d = 0;
        malloc_Int_Vector(&arc_d, arc_n);

        for (int j = 0; j < arc_n; j++){
            int n_epoch = arc[j][1] - arc[j][0];
            if (n_epoch<10){
                for (int k = arc[j][0]; k < arc[j][1]+1; k++){
                    obs_temp->P1[k][prn] = 0; obs_temp->P2[k][prn] = 0; obs_temp->L1[k][prn] = 0; obs_temp->L2[k][prn] = 0;
                    L6[k][prn] = 0; Li[k][prn] = 0; Nw[k][prn] = 0;
                }
                arc_d[num_of_arc_d] = j;
                num_of_arc_d++;
            }
        }

        //删除掉arc_d中的行
        for (int j = num_of_arc_d-1 ; j >= 0; j--){
            deleteRowAndReplace(&arc, &arc_n, arc_d[j]);
        }

        //----mw detect cycle slip------------------
        int j = 0;

        while (j<arc_n){
            //----first epoch check----------
            int e = arc[j][0];
            while (1){
                if (e+1 == arc[j][1] || e == arc[j][1]){
                    break;
                }
                double fir = Nw[e][prn]; double sec = Nw[e+1][prn]; double thi = Nw[e+2][prn];
                double firl = Li[e][prn]; double secl = Li[e+1][prn]; double thil = Li[e+2][prn];
                double sub = fabs(fir - sec); double sub2 = fabs(sec - thi);
                double subl = fabs(firl - secl); double subl2 = fabs(secl - thil);

                if (sub>1 || sub2>1 || subl>1 || subl2>1){
                    L6[e][prn] = 0; Li[e][prn] = 0; Nw[e][prn] = 0;
                    obs_temp->L1[e][prn] = 0; obs_temp->L2[e][prn] = 0; obs_temp->P1[e][prn] = 0; obs_temp->P2[e][prn] = 0;
                    e++;
                    arc[j][0] = e;
                }else{
                    arc[j][0] = e;
                    break;
                }
            }

            //----detect------------------
            if (arc[j][1] - arc[j][0]<10){
                for (int k = arc[j][0]; k < arc[j][1]+1; k++){
                    obs_temp->P1[k][prn] = 0; obs_temp->P2[k][prn] = 0; obs_temp->L1[k][prn] = 0; obs_temp->L2[k][prn] = 0;
                    L6[k][prn] = 0; Li[k][prn] = 0; Nw[k][prn] = 0;
                }
                //删除掉arc_final[j]这一行
                deleteRowAndReplace(&arc, &arc_n, j);
                continue;
            }

            double* ave_N; double* sigma; double* sigma2;
            malloc_Double_Vector(&ave_N, (arc[j][1] - arc[j][0] + 1));
            malloc_Double_Vector(&sigma, (arc[j][1] - arc[j][0] + 1));
            malloc_Double_Vector(&sigma2, (arc[j][1] - arc[j][0] + 1));

            ave_N[0] = Nw[arc[j][0]][prn];
            sigma[0] = 0; sigma2[0] = 0;
            int count = 1;
            double T = 0; double I1 = 0; double I2 = 0;

            //----------------------check epoch k+1
            for (int k = arc[j][0]+1; k < arc[j][1]; k++){
                ave_N[count] = ave_N[count-1] + (Nw[k][prn] - ave_N[count-1])/(count+1);
                sigma2[count] = sigma2[count-1] + ((Nw[k][prn] - ave_N[count-1])*(Nw[k][prn] - ave_N[count-1]) - sigma2[count-1])/(count+1);
                sigma[count] = sqrt(sigma2[count]/(count+1));
                T = fabs(Nw[k+1][prn] - ave_N[count]);
                I1 = fabs(Li[k+1][prn] - Li[k][prn]);

                //-------------------------no cycle slip
                if (T<4*sigma[count]&&I1<0.28){
                    count++;
                    continue;
                }else{
                    //---------------------arc end
                    if (k+1 == arc[j][1]){
                        if(k+1-arc[j][0]>10){
                            L6[k+1][prn] = 0; Li[k][prn] = 0; Nw[k+1][prn] = 0;
                            obs_temp->L1[k+1][prn] = 0; obs_temp->L2[k+1][prn] = 0; obs_temp->P1[k+1][prn] = 0; obs_temp->P2[k+1][prn] = 0;
                            arc[j][1] = k;
                        }else{//------delete scatter epoches
                            for (int l = arc[j][0]; l<k+1+1;l++){
                                L6[l][prn] = 0; Li[l][prn] = 0; Nw[l][prn] = 0;
                                obs_temp->L1[l][prn] = 0; obs_temp->L2[l][prn] = 0; obs_temp->P1[l][prn] = 0; obs_temp->P2[l][prn] = 0;
                            }
                            //删去第j行
                            deleteRowAndReplace(&arc, &arc_n, j);//arc_final_n--已经在函数中操作了
                            j--;
                        }
                        break;
                    }

                    I2 = fabs(Li[k+2][prn] - Li[k+1][prn]);

                    //还没debug--这个函数剩余部分
                    if (fabs(Nw[k+2][prn] - Nw[k+1][prn])<1&&I2<1){ //-----------cycle slip
                        if (k+1-arc[j][0]>10) {
                            //在arc中插入一个k。具体做法：第1部分-arc的前j-1行；第2部分-arc的第j行第一个元素和k；第3部分-k+1到arc第j行的第2个元素；第4部分-arc的第j+1行到最后一行
                            InsertRowInArc(&arc, &arc_n, j, k,1);
                        }else{
                            for (int l = arc[j][0]; l < k+1; l++){
                                L6[l][prn] = 0; Li[l][prn] = 0; Nw[l][prn] = 0;
                                obs_temp->L1[l][prn] = 0; obs_temp->L2[l][prn] = 0; obs_temp->P1[l][prn] = 0; obs_temp->P2[l][prn] = 0;
                            }
                            arc[j][0] = k+1;
                            j--;
                        }
                    }else{  //-----------------gross error
                        if (k+1-arc[j][0]>10){
                            L6[k+1][prn] = 0; Li[k][prn] = 0; Nw[k+1][prn] = 0;
                            obs_temp->L1[k+1][prn] = 0; obs_temp->L2[k+1][prn] = 0; obs_temp->P1[k+1][prn] = 0; obs_temp->P2[k+1][prn] = 0;
                            InsertRowInArc(&arc, &arc_n, j, k,2);
                        }else{
                            for (int l = arc[j][0]; l<k+2;l++){
                                L6[l][prn] = 0; Li[l][prn] = 0; Nw[l][prn] = 0;
                                obs_temp->L1[l][prn] = 0; obs_temp->L2[l][prn] = 0; obs_temp->P1[l][prn] = 0; obs_temp->P2[l][prn] = 0;
                            }
                            arc[j][0] = k+2;
                            j--;
                        }
                    }
                    break;
                }
            }
            j++;
        }
        //P4(:,i)=obs.P1(:,i)-obs.P2(:,i);
        for (int j = 0; j < 2880; j++){
            P4[j][prn] = obs_temp->P1[j][prn] - obs_temp->P2[j][prn];
        }

        //L4(:,i)=(c/f1)*obs.L1(:,i)-(c/f2)*obs.L2(:,i);
        for (int j = 0; j < 2880; j++){
            L4[j][prn] = (c/f1)*obs_temp->L1[j][prn] - (c/f2)*obs_temp->L2[j][prn];
        }

        //--------smoothing-------------------------
        for (int j = 0; j<arc_n; j++){
            int t = 2;
            for (int k = arc[j][0]+1; k<arc[j][1]+1;k++){
                P4[k][prn] = P4[k][prn]/t + (P4[k-1][prn]+L4[k-1][prn]-L4[k][prn])*(t-1)/t;
                t++;
            }
            //P4(arc(j,1):arc(j,1)+4,i)=0;
            for (int k = arc[j][0]; k<arc[j][0]+5;k++){
                P4[k][prn] = 0;
            }
        }

        //--------remove bad P4---------------------
        arc = Get_arc(P4, prn, &arc_n);

        for (int j = 0; j<arc_n; j++){
            double ave = mean(P4, arc[j][0], arc[j][1], prn);

            if (fabs(ave)>10){
                for (int kk = arc[j][0]; kk<arc[j][1]+1;kk++){
                    P4[kk][prn] = 0;
                }
            }
        }
    }

    //打印P4矩阵
//    for (int i = 0; i < 2880; i++){
//        printf(" 第%d行： ", i+1);
//        for (int j = 0; j < 32; j++){
//            printf(" %f ", P4[i][j]);
//        }
//        printf("\n");
//    }

    return P4;
}

int** Get_arc(double** L6, int prn, int* parc_n){
    int* arc; int index_of_arc = 0;
    malloc_Int_Vector(&arc, 2880);

    //arc赋初值为0
    for (int i = 0; i < 2880; i++){
        arc[i] = 0;
    }

    for (int i = 0; i < 2880; i++){
        if (i == 2880 - 1){
            if (L6[i][prn] != 0){
                arc[index_of_arc] = i;
                index_of_arc++;
            }
            continue;
        }
        if (i == 0 && L6[i][prn] != 0){
            arc[index_of_arc] = i;
            index_of_arc++;
        }
        //if array(i)==0&&array(i+1)~=0
        if (L6[i][prn] == 0 && L6[i+1][prn] != 0){
            arc[index_of_arc] = i+1;
            index_of_arc++;
            continue;
        }
        //if array(i)~=0&&array(i+1)==0
        if (L6[i][prn] != 0 && L6[i+1][prn] == 0){
            arc[index_of_arc] = i;
            index_of_arc++;
            continue;
        }
    }

    int returnSize;
    int** result = ConvertVectorTo2DMatrix(arc, index_of_arc, &returnSize);

    *parc_n = returnSize;//返回arc的行数

    return result;
}

// 功能：删除二维数组的某行，并返回删除后的数组
// 入参：1. arc：指向二维数组的指针 2. n：指向数组行数的指针 3. row：要删除的行的索引
void deleteRowAndReplace(int ***arc, int *n, int row) {
    int newN = *n - 1; // 新数组的行数
    int **tempArc;
    malloc_Int_2D(&tempArc, newN, 2);

    // 复制除了要删除的行之外的所有行到临时数组
    for (int i = 0, j = 0; i < *n; i++) {//这里应该从后面往前面删，不然会删一个，序号变一个
        if (i != row) {
            tempArc[j][0] = (*arc)[i][0];
            tempArc[j][1] = (*arc)[i][1];
            j++;
        }
    }

    // 释放原数组的内存
    for (int i = 0; i < *n; i++) {
        free((*arc)[i]);
    }
    free(*arc);

    // 将原数组指针指向新数组
    *arc = tempArc;
    *n = newN; // 更新外部n的值
}

double mean(double** array, int start, int end, int prn){
    double sum = 0;
    for (int i = start; i < end; i++){
        sum += array[i][prn];
    }
    return sum / (end - start);
}

void save_P4_To_Bin_File(double **P4, int dim1, int dim2, const char *filename) {
    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
        perror("Failed to open file");
        exit(1);
    }

    for (int i = 0; i < dim1; ++i) {
        fwrite(P4[i], sizeof(double), dim2, file);
    }

    fclose(file);
}

void load_P4_From_Bin_File(double ***P4, int dim1, int dim2, const char *filename) {
    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
        perror("Failed to open file");
        exit(1);
    }

    *P4 = (double **)malloc(dim1 * sizeof(double *));
    for (int i = 0; i < dim1; ++i) {
        (*P4)[i] = (double *)malloc(dim2 * sizeof(double));
        fread((*P4)[i], sizeof(double), dim2, file);
    }

    fclose(file);
}

void DCB_Estimation(){
    printf("----------Step 5: starting DCB Estimation!----------\n");
    //1. 先用一个函数得到一个去重的doy_diff数组
    //2. 然后遍历这个doy_diff数组，对每一个doy，load对应的SP3文件（卫星坐标）
    //3. 将取出的sate_xyz, sitesInfo, doy, SDCB_REF, 和order作为入参传入get_MDCB。返回量包括DCB_R, DCB_S, IONC三个估计量结果数组。
    int* doy_diff; int doy_diff_size;
    removeSameElementAndReturn_Int(sitesInfo.doy, &doy_diff, r_file_num, &doy_diff_size);

    for (int i = 0; i < doy_diff_size; i++){
        int doy = doy_diff[i];
        double*** sate_xyz;
        malloc_double_3D(&sate_xyz, 3, 2880, N_SAT);

        char* load_sp3_file_path = get_load_sp3_pathname(doy);
        load_Sp3_From_Bin_File(&sate_xyz, 3, 2880, N_SAT, load_sp3_file_path);
        double* DCB_R; double* DCB_S; double* IONC;
        int DCB_R_size = 0; int DCB_S_size = 0; int IONC_size = 0;
        get_DCB(&DCB_R, &DCB_R_size, &DCB_S, &DCB_S_size, &IONC, &IONC_size,doy, sate_xyz, sitesInfo, sDCB_REF, order);



        char* DCB_R_final_result_filename = generate_final_result_file_pathname(doy, "DCB_R", m_result_path);
        char* DCB_S_final_result_filename = generate_final_result_file_pathname(doy, "DCB_S", m_result_path);
        char* IONC_final_result_filename = generate_final_result_file_pathname(doy, "IONC", m_result_path);
        write_vector_to_file(DCB_R_final_result_filename, DCB_R, DCB_R_size);
        write_vector_to_file(DCB_S_final_result_filename, DCB_S, DCB_S_size);
        write_vector_to_file(IONC_final_result_filename, IONC, IONC_size);
//        printf("Successfully saved %s to file\n", DCB_R_final_result_filename);
//        printf("Successfully saved %s to file\n", DCB_S_final_result_filename);
//        printf("Successfully saved %s to file\n", IONC_final_result_filename);
    }
    printf("----------Step 5: completings !----------\n");
}

// 返回一个包含文件夹中文件的名称的字符串数组，注意入参是const char
// 使用前请先使用countFilesInDirectory函数获取文件数量
// 使用前请提前声明一个对应尺寸的字符串数组
char** getFileNamesInDirectory(const char *path, char **fileNames, int num) {
    DIR *dir;
    struct dirent *ent;
    int i = 0;
    if ((dir = opendir(path)) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            if (ent->d_name[0] == '.'){
                continue;
            }
            strcpy(fileNames[i], ent->d_name);
            i++;
        }
        closedir(dir);
    } else {
        perror("Failed to open directory");
    }
    return fileNames;
}

// 函数声明，用于创建新的二维数组
// 使用前请先计算arr的行数，例如：int size = sizeof(arr) / sizeof(arr[0]);
// 使用前请先声明一个对应尺寸的二维数组，但不需要赋值（因为函数内最后会通过指针修改这个值），例如： int returnSize;
int** ConvertVectorTo2DMatrix(int *arr, int size, int *returnSize) {
    // 动态分配最大可能需要的空间，即size/2个子数组
    int** result = (int**)malloc((size / 2) * sizeof(int*));
    int count = 0; // 用于记录最终二维数组的行数

    for (int i = 0; i < size; i += 2) {
        // 跳过所有两个元素都是0的子数组
        if (i + 1 < size && !(arr[i] == 0 && arr[i + 1] == 0)) {
            result[count] = (int*)malloc(2 * sizeof(int)); // 为子数组分配空间
            result[count][0] = arr[i];   // 存储当前元素和下一个元素
            result[count][1] = arr[i + 1];
            count++; // 增加二维数组的行数
        }
    }

    *returnSize = count; // 通过指针参数返回二维数组的实际行数
    return result; // 返回创建的二维数组
}

//
void InsertRowInArc(int ***arc, int *n, int row, int k, int k_shift) {
    int newN = *n + 1; // 新数组的行数
    int **tempArc = (int **)malloc(newN * sizeof(int *));
    for (int i = 0; i < newN; i++) {
        tempArc[i] = (int *)malloc(2 * sizeof(int));
    }

    //第1部分：arc的前j-1行
    for (int i = 0; i < row; i++) {
        tempArc[i][0] = (*arc)[i][0];
        tempArc[i][1] = (*arc)[i][1];
    }

    //第2部分：arc的第j行第一个元素和k
    tempArc[row][0] = (*arc)[row][0];
    tempArc[row][1] = k;

    //第3部分：k+1到arc第j行的第2个元素
    tempArc[row+1][0] = k+k_shift;
    tempArc[row+1][1] = (*arc)[row][1];

    //第4部分：arc的第j+1行到最后一行
    for (int i = row+2; i < newN; i++) {
        tempArc[i][0] = (*arc)[i-1][0];
        tempArc[i][1] = (*arc)[i-1][1];
    }

    // 释放原数组的内存
    for (int i = 0; i < *n; i++) {
        free((*arc)[i]);
    }
    free(*arc);

    // 将原数组指针指向新数组
    *arc = tempArc;
    *n = newN; // 更新外部n的值
}

void removeSameElementAndReturn_Char(char** vector, char*** new_vector, int vector_size, int* new_size){
    //给doys申请n个int类型的内存，并赋初值为0
    int new_size_temp = 0;

    //先申请一个跟原数组一样尺寸的数组，用于存储去重后的数组
    char** new_vector_temp;
    malloc_Char_Vector(&new_vector_temp, vector_size, 64);

    //得到去重后的数组。但此时还不能直接返回，因为去重后的数组中可能有空字符串
    strcpy(new_vector_temp[new_size_temp], vector[0]);
    new_size_temp++;
    for (int i = 0; i < vector_size; i++) {
        //先判断new_vector_temp中是否已经有了vector[i]，如果有了就跳过
        int flag = 0;
        for (int j = 0; j < new_size_temp; j++) {
            if (strcmp(vector[i], new_vector_temp[j]) == 0) {
                flag = 1;
                break;
            }
        }
        //如果new_vector_temp中没有vector[i]，就将vector[i]加入new_vector_temp
        if (flag == 0) {
            strcpy(new_vector_temp[new_size_temp], vector[i]);
            new_size_temp++;
        }
    }

    //再申请一个新的数组，用于存储去掉空字符串后的数组
    char** new_vector_temp_2;
    malloc_Char_Vector(&new_vector_temp_2, new_size_temp, 64);
    for (int i = 0; i < new_size_temp; i++){
        strcpy(new_vector_temp_2[i], new_vector_temp[i]);
    }

    //将new_vector的指针指向new_vector_temp_2
    *new_vector = new_vector_temp_2;
    //修改外部new_size的值
    *new_size = new_size_temp;

    free(new_vector_temp);
}

void removeSameElementAndReturn_Int(int* vector, int** new_vector, int vector_size, int* new_size){
    // 申请一个临时数组用于存储去重后的数组
    int* new_vector_temp = (int*)malloc(vector_size * sizeof(int));
    int new_size_temp = 0;

    // 得到去重后的数组，但此时还不能直接返回，因为去重后的数组中可能有重复的元素
    new_vector_temp[new_size_temp] = vector[0];
    new_size_temp++;
    for (int i = 1; i < vector_size; i++) {
        // 先判断new_vector_temp中是否已经有了vector[i]，如果有了就跳过
        int flag = 0;
        for (int j = 0; j < new_size_temp; j++) {
            if (vector[i] == new_vector_temp[j]) {
                flag = 1;
                break;
            }
        }
        // 如果new_vector_temp中没有vector[i]，就将vector[i]加入new_vector_temp
        if (flag == 0) {
            new_vector_temp[new_size_temp] = vector[i];
            new_size_temp++;
        }
    }

    // 申请一个新的数组，用于存储去掉重复元素后的数组
    int* new_vector_temp_2 = (int*)malloc(new_size_temp * sizeof(int));
    memcpy(new_vector_temp_2, new_vector_temp, new_size_temp * sizeof(int));

    // 将new_vector的指针指向new_vector_temp_2
    *new_vector = new_vector_temp_2;
    // 修改外部new_size的值
    *new_size = new_size_temp;

    // 释放内存
    free(new_vector_temp);
}

char* get_load_sp3_pathname(int doy){
    char* filepathname = (char *) malloc(64 * sizeof(char));

    //将s_opath赋给filepathname
    strcpy(filepathname, s_opath);

#ifdef LINUX
    strcat(filepathname, "/");
#elif defined(WINDOWS)
    //加入"\\"
    strcat(filepathname, "\\");
#else
#endif

    //加入"20"
    strcat(filepathname, "20");

    char doy_str[3];
    snprintf(doy_str, sizeof(doy_str), "%d", doy);
    strncat(filepathname, doy_str, 2);

    //在字符串中添加"_"
    strcat(filepathname, "_");
    // 提取doy的后三位数字
    int last_three_digits = doy % 1000;
    // 转换后三位数字为字符串
    char doy_str_2[4]; // 考虑到最大值为999，因此最多需要4个字符（包括终止符'\0'）
    sprintf(doy_str_2, "%d", last_three_digits);
    strcat(filepathname, doy_str_2);

    //加入"sp3.dat"
    strcat(filepathname, "sp3.dat");

    return filepathname;
}

void get_DCB(double** DCB_R, int* DCB_R_size, double** DCB_S, int* DCB_S_size, double** IONC, int* IONC_size, int doy, double*** sate_xyz, SitesInfo sitesInfo, SDCB_REF sDCB_REF, int order) {
    char **stations;
    malloc_Char_Vector(&stations, r_file_num, 4);

    for (int i = 0; i < r_file_num; i++) {
        strncpy(stations[i], &sitesInfo.name[i][0], 4);
        stations[i][4] = '\0';
    }

    double **x = sate_xyz[0];
    double **y = sate_xyz[1];
    double **z = sate_xyz[2];

    int n_r = 0; int n_s = N_SAT;
    char **P4_load_filepath = get_doy_P4_filepath(doy, &n_r);
    sort_filenames(P4_load_filepath, n_r);

    //--check the number of each satellite's observations
    int *PRN;
    malloc_Int_Vector(&PRN, N_SAT);

    //为DCB_S申请内存
    *DCB_S = (double *) malloc(N_SAT * sizeof(double));

    //计算每个PRN的所有P4值
    for (int i = 0; i < n_r; i++) {
        double **P4;
        //为P4申请内存并赋初值为0
        malloc_Double_2D(&P4, 2880, N_SAT);

        load_P4_From_Bin_File(&P4, 2880, N_SAT, P4_load_filepath[i]);

        for (int j = 0; j < N_SAT; j++) {
            for (int k = 0; k < 2880; k++) {
                if (P4[k][j] != 0) {
                    PRN[j]++;
                }
            }
        }
    }

    // 构建最小二乘估计的矩阵
    //printf("Starting to build the matrix for the least squares estimation...\n");
    int est_num = (order+1)*(order+1)*12+n_s+n_r;

    double** B = NULL; double* l = NULL;

    double* C = (double *) malloc(est_num * sizeof(double));
    //赋初值为0
    for (int i = 0; i < est_num; i++){
        if (i>n_r-1&&i<n_r+n_s){
            C[i] = 1;
        }else{
            C[i] = 0;
        }
    }

    double Wx = 0; int B_row_number = 0; int l_row_number = 0;

    for (int i=0; i<n_r;i++){
        double **P4;
        //为P4申请内存并赋初值为0
        malloc_Double_2D(&P4, 2880, N_SAT);
        P4 = (double **) malloc(2880 * sizeof(double *));
        for (int j = 0; j < 2880; j++) {
            P4[j] = (double *) malloc(N_SAT * sizeof(double));
            for (int k = 0; k < N_SAT; k++) {
                P4[j][k] = 0;
            }
        }

        load_P4_From_Bin_File(&P4, 2880, N_SAT, P4_load_filepath[i]);
        //取P4_load_filepath[i]的前4个字符，赋给site
        char site[5];
        strncpy(site, P4_load_filepath[i], 4);
        site[4] = '\0';

        //找到stations里面同时满足site和doy的索引
        int index = 0;
//        for (int j = 0; j < r_file_num; j++){
//            if (strcmp(stations[j], site) == 0&&sitesInfo.doy[j] == doy){
//                index = j;
//                break;
//            }
//        }

        double sx = sitesInfo.coor[index][0];
        double sy = sitesInfo.coor[index][1];
        double sz = sitesInfo.coor[index][2];

        int ith = i;

        double** sN; double* sL;

        int sN_row_number = 0; int sL_number = 0;
        get_Matrix(P4, x, y, z, sx, sy, sz, n_s, n_r, ith, order, &sN, &sL, est_num, &sN_row_number, &sL_number);

        //把sN(n*337二维)的结果纵向添加进B矩阵(m*337二维)，得到((m+n)*337)的矩阵
        addMartixDownToAnother(&B, &sN, est_num, &B_row_number, &l_row_number, sN_row_number);
        addVectorToVector(&l, &sL, &l_row_number, sL_number);
    }


    // 添加零均值条件
    addRowTo2DMatrixNewLine(&B, &C, &B_row_number, est_num);
    appendDoubleToVector(&l, &l_row_number, Wx);

    // 求解最小二乘R = A\B
    double* R = solve_least_squares(B, l, B_row_number, est_num, l_row_number);

    double* DCB_R_temp = getDCB_R(R, n_r);
    double* DCB_S_temp = getDCB_S(R, n_r, n_s);
    double* IONC_temp = getIONC(R, n_r, n_s);

    *DCB_R = DCB_R_temp;
    *DCB_R_size = n_r;
    *DCB_S = DCB_S_temp;
    *DCB_S_size = n_s;
    *IONC = IONC_temp;
    *IONC_size = (order+1)*(order+1)*12;
}

//找到对应doy的各个站所有P4文件
char** get_doy_P4_filepath(int doy, int* n_r){
    //遍历m_p4_path文件夹中的文件名，取出每个文件名的第5到第9位，转换成int，如果与doy相同，就将这个文件名存入一个字符串数组中
    int p4_num = countFilesInDirectory(m_p4_path);
    char** p4_file_names;
    malloc_Char_Vector(&p4_file_names, p4_num, 64);
    p4_file_names = getFileNamesInDirectory(m_p4_path, p4_file_names, p4_num);

    char** doy_p4_file_names;
    //分配足够内存
    malloc_Char_Vector(&doy_p4_file_names, p4_num, 64);

    int doy_p4_num = 0;
    for(int i = 0; i < p4_num; i++){
        char* doy_str = (char *) malloc(64 * sizeof(char));
        strncpy(doy_str, p4_file_names[i]+4, 5);
        doy_str[5] = '\0';
        int doy_temp = atoi(doy_str);
        if (doy_temp == doy){
            //将这个文件名存入一个字符串数组中
            strcpy(doy_p4_file_names[doy_p4_num], p4_file_names[i]);
            doy_p4_num++;
        }
    }

    char** doy_p4_file_names_return;
    malloc_Char_Vector(&doy_p4_file_names_return, doy_p4_num, 64);

    //把doy_p4_file_names中的前doy_p4_num个元素复制到doy_p4_file_names_return中
    for (int i = 0; i < doy_p4_num; i++){
        strcpy(doy_p4_file_names_return[i], doy_p4_file_names[i]);
    }

    //释放doy_p4_file_names的内存
    for (int i = 0; i < doy_p4_num; i++){
        free(doy_p4_file_names[i]);
    }
    free(doy_p4_file_names);

    *n_r = doy_p4_num;

    //return
    return doy_p4_file_names_return;
}

void get_Matrix(double** P4, double** x, double** y, double**z, double sx, double sy, double sz, int n_s, int n_r, int ith, int order, double*** sN, double** sL, int est_num, int* M_row_number, int* l_size_number){
    double** M = NULL ; double* l = NULL; int l_size = 0;
    //给M分配1行337列的内存
    int M_row_num = 0; int M_col_num = est_num;
//    M = (double **) malloc(M_row_num * sizeof(double *));
//    for (int i = 0; i < M_row_num; i++) {
//        M[i] = (double *) malloc(M_col_num * sizeof(double));
//        for (int j = 0; j < 337; j++) {
//            M[i][j] = 0.0; // 初始化为零
//        }
//    }

    double* BL = XYZ2BLH(sx, sy, sz);
    double sb = BL[0]; double sl = BL[1];

    for (int i = 0 ;i<(24/period); i++){
        for (int j = 0; j<N_SAT; j++){
            for (int k=240*i; k<240*(i+1);k++) {
                if (P4[k][j]==0){
                    continue;
                }

                int est_num = (order+1)*(order+1)*12+n_s+n_r;

                double* M_col;
                malloc_Double_Vector(&M_col, est_num);

                double* EA = get_EA(sx, sy, sz, x[k][j]*1000, y[k][j]*1000, z[k][j]*1000);
                double E = EA[0]; double A = EA[1];

                double IPPz = asin(6378137*sin(0.9782*(pi/2-E))/(6378137+506700));
                double t_r = 30*k*pi/43200;

                double* BS = get_IPP(E, A, sb, sl, IPPz, t_r);
                double b = BS[0]; double s = BS[1];

                M_col[ith] = (-9.52437)*cos(IPPz);
                M_col[n_r+j] = (-9.52437)*cos(IPPz);
                int st = (order+1)*(order+1)*i+n_r+N_SAT;
                int ed = (order+1)*(order+1)*(i+1)+n_r+N_SAT;


                get_Coef(&M_col, st, ed, b, s, order);


                //把M_col添加到二维数组M的下一行
                addRowTo2DMatrixNewLine(&M, &M_col, &M_row_num, M_col_num);

                //把P4[k][j]*(-9.52437)*cos(IPPz)加到l中
                appendDoubleToVector(&l, &l_size, P4[k][j]*(-9.52437)*cos(IPPz));
            }
        }
    }
    //返回M的行数和l的个数
    *M_row_number = M_row_num;
    *l_size_number = l_size;
    //返回M和l
    *sN = M;
    *sL = l;
}


double* get_IPP(double E, double A, double sb, double sl, double IPPz, double t_r){
    double* BS = (double *) malloc(2 * sizeof(double));
    BS[0] = 0; BS[1] = 0;

    double t = pi/2 - E - IPPz;
    BS[0] = asin(sin(sb)*cos(t)+cos(sb)*sin(t)*cos(A));
    BS[1] = sl + asin(sin(t)*sin(A)/cos(BS[0]));
    BS[1] = BS[1] + t_r - pi;

    return BS;
}

void get_Coef(double** M_col, int st, int ed, double b, double s, int order){
    double* cof_P = (double *) malloc((order+1)*(order+1) * sizeof(double));

    double* ms;
    ms = (double *) malloc(order * sizeof(double));
    for (int i = 0; i < order; i++){
        ms[i] = s*(i+1);
    }

    int i = 0; double N = 0;
    double x = sin(b);
    for (int n = 0; n < order+1; n++){
        double* P = get_legendre(n, x);
        for (int m = 0; m<n+1; m++){
            if (m == 0){
                N = norm(n,m);
                cof_P[i] = P[m]*N;
            }else{
                N = norm(n,m);
                cof_P[i] = P[m]*N*cos(ms[m-1]);
                i++;
                N = norm(n,m);
                cof_P[i] = P[m]*N*sin(ms[m-1]);
            }
            i++;
        }
    }

    //把cof_P添加到M_col的st到ed的位置
    for (int i = st; i<ed; i++){
        (*M_col)[i] = cof_P[i-st];
    }

}

double* get_legendre(int n, double x){
    double* P = (double *) malloc((n+1) * sizeof(double));
    for (int i = 0; i < n+1; i++){
        //如果i为奇数
        if (i%2 == 1) {
            P[i] = -(std::assoc_legendre(n, i, x));
        }else{
            P[i] = std::assoc_legendre(n, i, x);
        }
    }

    return P;
}

// 函数说明：自己照着Matlab的lagendre函数写的，用不了
double* legendre(int n, double x){
    double* P; int P_size = 0;
    malloc_Double_Vector(&P, n+1);

    if (n==0){
        double y = 1;
        P[0] = y;
        P_size++;
    }

    double* rootn;
    rootn = (double *) malloc((2*n+1) * sizeof(double));
    for (int i = 0; i < 2*n+1; i++){
        rootn[i] = sqrt(i);
    }

    double s = sqrt(1-x*x);

    double twocot = (-2)*x/s;

    double sn = pow((-s),n);

    double tol = sqrt(DBL_MIN);

    if (s>0&&fabs(sn)<=tol){

    }

    if (x!=1&&fabs(sn)>=tol){
        double* d = (double *) malloc(n * sizeof(double));
        for (int i = 0; i < n; i++){
            d[i] = 2+2*i;
        }

        double product = 1;
        for (int i = 0; i < n; i++){
            product = product*(1-(1/d[i]));
        }

        P[n] = sqrt(product)*sn;
        P_size++;
        P[n-1] = P[n]*twocot*n/rootn[2*n];
        P_size++;

        for (int m=n-2; m>=0; m--){
            P[m] = (P[m+1]*twocot*(m+1)-P[m+2]*rootn[n+m+2]*rootn[n-m-1]/(rootn[n+m+1]*rootn[n-m]));
            P_size++;
        }
    }

    if (s==0){
    }

    for (int m = 0; m < n-1; m++){
        P[m+1] = P[m+1]*rootn[2*m];//这行不确定有没有问题
    }
    P[n] = P[n]*rootn[2*n];

    return P;
}

double norm(int n, int m){
    double N = 0;
    if (m==0){
        N = sqrt(factorial(n-m)*(2*n+1)/factorial(n+m));
    }else{
        N = sqrt(factorial(n-m)*(4*n+2)/factorial(n+m));
    }
    return N;
}

double factorial(int n) {
    if (n < 0) {
        // 阶乘未定义于负数
        return 0;
    }
    int result = 1;
    for (int i = 2; i <= n; ++i) {
        result *= i;
    }
    return result;
}

// 添加新行到M
void addRowTo2DMatrixNewLine(double ***M, double **M_sol, int *numRows, int numCols) {
    // 增加行数
    (*numRows)++;

    // 重新分配内存
    *M = (double **)realloc(*M, (*numRows) * sizeof(double *));
    if (*M == NULL) {
        printf("内存分配失败\n");
        exit(1);
    }

    // 分配每行的内存
    (*M)[*numRows - 1] = (double *)malloc(numCols * sizeof(double));
    if ((*M)[*numRows - 1] == NULL) {
        printf("内存分配失败\n");
        exit(1);
    }

    // 将新的一维数组添加到M中
    (*M)[*numRows - 1] = *M_sol;
}

// 向一维数组尾部添加一个double元素
void appendDoubleToVector(double **arr, int *size, double value) {
    // 增加数组大小
    (*size)++;

    // 重新分配内存
    *arr = (double *)realloc(*arr, (*size) * sizeof(double));
    if (*arr == NULL) {
        printf("内存分配失败\n");
        exit(1);
    }

    // 将值添加到数组尾部
    (*arr)[(*size) - 1] = value;
}

void addMartixDownToAnother(double ***B, double*** sN, int est_num, int* B_row_number, int* l_row_number, int sN_row_number){
    // 增加行数
    (*B_row_number) += sN_row_number;

    // 重新分配内存
    *B = (double **)realloc(*B, (*B_row_number) * sizeof(double *));
    if (*B == NULL) {
        printf("内存分配失败\n");
        exit(1);
    }

    // 分配每行的内存
    for (int i = 0; i < sN_row_number; i++){
        (*B)[*B_row_number - sN_row_number + i] = (double *)malloc(est_num * sizeof(double));
        if ((*B)[*B_row_number - sN_row_number + i] == NULL) {
            printf("内存分配失败\n");
            exit(1);
        }
    }

    // 将新的二维数组添加到B中
    for (int i = 0; i < sN_row_number; i++){
        (*B)[*B_row_number - sN_row_number + i] = (*sN)[i];
    }
}

void addVectorToVector(double** l, double** sL, int* l_row_number, int sL_size){
    // 增加行数
    (*l_row_number) += sL_size;

    // 重新分配内存
    *l = (double *)realloc(*l, (*l_row_number) * sizeof(double));
    if (*l == NULL) {
        printf("内存分配失败\n");
        exit(1);
    }

    // 将新的一维数组添加到l中
    for (int i = 0; i < sL_size; i++){
        (*l)[(*l_row_number) - sL_size + i] = (*sL)[i];
    }
}

double* getDCB_R(double* R, int n_r){
    double* DCB_R = (double *) malloc(n_r * sizeof(double));
    for (int i = 0; i < n_r; i++){
        DCB_R[i] = R[i]*10e8/c;
    }
    return DCB_R;
}

double* getDCB_S(double* R, int n_r, int n_s){
    double* DCB_S = (double *) malloc(n_s * sizeof(double));
    for (int i = 0; i < n_s; i++){
        DCB_S[i] = R[i+n_r]*10e8/c;
    }
    return DCB_S;
}

double* getIONC(double* R, int n_r, int n_s){
    double* IONC = (double *) malloc(((order+1)*(order+1)*12) * sizeof(double));
    for (int i = 0; i < (order+1)*(order+1)*12; i++){
        IONC[i] = R[i+n_s+n_r];
    }
    return IONC;
}

void write_vector_to_file(char *filename, double* vector, int size) {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Error opening file.\n");
        return;
    }
    for (int i = 0; i < size; ++i) {
        fprintf(fp, "%lf\n", vector[i]);
    }
    fclose(fp);
}

char* generate_final_result_file_pathname(int doy, const char* type, const char* m_result_path){
    char* pathname = (char *) malloc(128 * sizeof(char));
    //添加m_result_path
    strcpy(pathname, m_result_path);
#ifdef LINUX
    strcat(pathname, "/");
#elif defined(WINDOWS)
    //添加"\\"
    strcat(pathname, "\\");
#else
#endif
    //添加type
    strcat(pathname, type);
    //添加doy
    char doy_str[5];
    sprintf(doy_str, "%d", doy);
    strcat(pathname, doy_str);
    //添加".txt"
    strcat(pathname, ".txt");

    return pathname;
}

// 为二维数组申请内存，并赋初值为0
void malloc_Double_2D(double*** array, int row, int col){
    *array = (double **) malloc(row * sizeof(double *));
    for (int i = 0; i < row; i++) {
        (*array)[i] = (double *) malloc(col * sizeof(double));
        for (int j = 0; j < col; j++) {
            (*array)[i][j] = 0.0; // 初始化为零
        }
    }
}

void malloc_Double_Vector(double** vector, int size){
    *vector = (double *) malloc(size * sizeof(double));
    for (int i = 0; i < size; i++){
        (*vector)[i] = 0.0;
    }
}

void malloc_Char_Vector(char*** vector, int size, int string_size){
    *vector = (char **) malloc(size * sizeof(char *));
    for (int i = 0; i < size; i++){
        (*vector)[i] = (char *) malloc(string_size * sizeof(char));
    }
}

void malloc_Int_Vector(int** vector, int size){
    *vector = (int *) malloc(size * sizeof(int));
    for (int i = 0; i < size; i++){
        (*vector)[i] = 0;
    }
}

void malloc_Int_2D(int*** array, int row, int col){
    *array = (int **) malloc(row * sizeof(int *));
    for (int i = 0; i < row; i++) {
        (*array)[i] = (int *) malloc(col * sizeof(int));
        for (int j = 0; j < col; j++) {
            (*array)[i][j] = 0; // 初始化为零
        }
    }
}

void malloc_double_3D(double**** array, int dim1, int dim2, int dim3){
    *array = (double ***) malloc(dim1 * sizeof(double **));
    for (int i = 0; i < dim1; i++) {
        (*array)[i] = (double **) malloc(dim2 * sizeof(double *));
        for (int j = 0; j < dim2; j++) {
            (*array)[i][j] = (double *) malloc(dim3 * sizeof(double));
            for (int k = 0; k < dim3; k++) {
                (*array)[i][j][k] = 0.0; // 初始化为零
            }
        }
    }
}

// 比较两个字符串，实现数字的自然排序，part1
int strnatcmp(const char *a, const char *b) {
    while (*a && *b) {
        if (isdigit(*a) && isdigit(*b)) {
            // 如果当前位置都是数字，则读取整个数字进行比较
            char *end_a;
            char *end_b;
            long num_a = strtol(a, &end_a, 10);
            long num_b = strtol(b, &end_b, 10);
            if (num_a < num_b) return -1;
            if (num_a > num_b) return 1;
            // 如果数字相等，继续比较数字后面的部分
            a = end_a;
            b = end_b;
        } else {
            // 普通字符的比较
            if (*a < *b) return -1;
            if (*a > *b) return 1;
            ++a;
            ++b;
        }
    }
    // 处理字符串长度不同的情况
    if (*a) return 1; // a 比 b 长
    if (*b) return -1; // b 比 a 长
    return 0;
}

// 用于qsort的比较函数包装（比较两个字符串，实现数字的自然排序，part2）
int compare(const void *a, const void *b) {
    const char *str_a = *(const char**)a;
    const char *str_b = *(const char**)b;
    return strnatcmp(str_a, str_b);
}

// 对字符串数组进行自然排序（比较两个字符串，实现数字的自然排序，part3）
void sort_filenames(char *files[], int count) {
    qsort(files, count, sizeof(char*), compare);
}

// 定义函数来转换二维数组B和一维数组l，并计算最小二乘解
double* solve_least_squares(double** B, double* l, int B_row_num, int B_col_num, int l_row_num) {
    cholmod_common Common, *cc;
    cholmod_sparse *B_chol;
    cholmod_dense *l_chol, *X_chol;
    double *B_chol_values_mat; int *B_chol_row_mat, *B_chol_col_mat;
    malloc_Double_Vector(&B_chol_values_mat, B_row_num * B_col_num);
    malloc_Int_Vector(&B_chol_row_mat, B_row_num * B_col_num);
    malloc_Int_Vector(&B_chol_col_mat, B_row_num * B_col_num);
    double *X, *X_values;
    int i, j;
    int nzz_max = 0;

    // 初始化CHOLMOD环境
    cc = &Common;
    cholmod_start(cc);

    // 将二维数组B转换为csc存储模式
    for (j = 0; j < B_col_num; j++) {
        for (i = 0; i < B_row_num; i++) {
            if (B[i][j] != 0) {
                // 如果B[i][j]非0，则添加到稀疏矩阵中
                B_chol_row_mat[nzz_max] = i;
                B_chol_col_mat[nzz_max] = j;
                B_chol_values_mat[nzz_max] = B[i][j];
                nzz_max++;
            }
        }
    }

    cholmod_triplet* triplet = cholmod_allocate_triplet(B_row_num,B_col_num,B_row_num*B_col_num,0,CHOLMOD_REAL,cc);
    int * triplet_i = (int *)(triplet->i);
    int * triplet_j = (int *)(triplet->j);
    double * triplet_x = (double *)(triplet->x);
    for (int ne=0; ne<nzz_max; ne++)
    {
        triplet_i[triplet->nnz] = B_chol_row_mat[ne];
        triplet_j[triplet->nnz] = B_chol_col_mat[ne];
        triplet_x[triplet->nnz] = B_chol_values_mat[ne];
        triplet->nnz++;
    }

    B_chol = cholmod_triplet_to_sparse(triplet, nzz_max, cc);
    cholmod_free_triplet(&triplet, cc);


    // 创建稀疏矩阵l_chol
    l_chol = cholmod_zeros(l_row_num, 1, CHOLMOD_REAL, cc);
    // 填充稀疏矩阵l_chol
    for (int ne=0; ne<l_row_num; ne++)
    {
        ((double *)(l_chol->x))[ne] = l[ne];
    }

    // 使用SuiteSparseQR计算最小二乘解
    X_chol = SuiteSparseQR_C_backslash_default(B_chol, l_chol, cc);

    // 从X_chol中提取解向量
    X_values = (double*)X_chol->x;
    malloc_Double_Vector(&X, l_row_num);
    for (i = 0; i < B_col_num; i++) {
        X[i] = X_values[i];
    }

    // 清理资源
    cholmod_free_sparse(&B_chol, cc);
    cholmod_free_dense(&l_chol, cc);
    cholmod_free_dense(&X_chol, cc);
    cholmod_finish(cc);

    return X; // 返回最小二乘解数组
}

// Linux创建文件
void CreateFolderIfNotExist_Linux(const char* folderPath) {
    struct stat st = {0};

    if (stat(folderPath, &st) == -1) {
        // 文件夹不存在，创建文件夹
        if (mkdir(folderPath, 0700) == 0) {
            printf("文件夹%s创建成功。\n",folderPath);
        } else {
            printf("文件夹%s创建失败。\n",folderPath);
        }
    }
}

#ifdef WINDOWS
void CreateFolderIfNotExist(const char* folderPath) {
    DWORD dwAttrib = GetFileAttributes(folderPath);

    if (dwAttrib == INVALID_FILE_ATTRIBUTES || !(dwAttrib & FILE_ATTRIBUTE_DIRECTORY)) {
        // 文件夹不存在，创建文件夹
        if (CreateDirectory(folderPath, NULL)) {
            printf("文件夹%s创建成功。\n",folderPath);
        } else {
            printf("文件夹%s创建失败。\n",folderPath);
        }
    }
}
#endif