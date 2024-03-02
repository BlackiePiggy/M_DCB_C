#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <math.h>

#define pi 3.14159265358979323846
#define c 299792458


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

// -----------------------------Setting--------------------------------------
const char* r_ipath = "D:\\projects\\M_DCB_C\\RINEX_files";// rinex
const char* r_opath = "D:\\projects\\M_DCB_C\\RINEX_output_files";//sprintf(sav_filename, "D:\\projects\\M_DCB\\RINEX_output_files\\observation_%s.csv", psitesInfo->name[i]);
const char* s_ipath = "D:\\projects\\M_DCB_C\\SP3_files";// sp3
const char* s_opath = "D:\\projects\\M_DCB_C\\SP3_output_files";//
const char* i_ipath = "D:\\projects\\M_DCB_C\\IONEX_files";// ionex
const char* i_opath = "D:\\projects\\M_DCB_C\\IONEX_output_files";// ionex
int lim = 10;// el
int order = 4;// order
int r_file_num;
double* DCB_rec;
double f1 = 1575.42e6;
double f2 = 1227.6e6;

//全局变量
SitesInfo sitesInfo;
Obs obs;
SDCB_REF sDCB_REF;

// -----------------------------Function--------------------------------------
int countFilesInDirectory(const char *folderPath);
//void free_usage(SitesInfo sitesInfo, int len);
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
int find_ionex_index(int doy);
int find_rinex_indexs(int* index, int doy);
void get_filenames(const char* folderPath, char** filenames, int file_num);
void read_single_ionex(const char* ionex_file, int i, char** sites_name, int doy_num, SDCB_REF sDCB_REF, double* DCB_rec);
void get_smoothed_P4(SitesInfo sitesInfo, double z_threshold, int flag);
void writeToFile(char* filename, Obs* obs, int num_sites);
void writeObsToFile(const char* filename, Obs* obs, int p1Rows, int p2Rows, int l1Rows, int l2Rows, int cols);
void readObsFromFile(const char* filename, Obs* obs);
void save_sp3_To_Bin_File(double ***sate_xyz, int dim1, int dim2, int dim3, const char *filename);
void load_Sp3_From_Bin_File(double ****sate_xyz, int dim1, int dim2, int dim3, const char *filename);
double* XYZ2BLH(double x, double y, double z);
double* get_EA(double sx, double sy, double sz, double x, double y, double z);
void cutobs(Obs *obs_temp, double*** sate_xyz, double sx, double sy, double sz, int lim);
void convertTo2D(int *arr, int size, int** result);
int** Get_arc(double** L6, int prn);
int** removeRows(int** arc, int n, int* arc_d, int d_size, int* new_size);
void deleteRowAndReplace(int ***arc, int *n, int row);

// -----------------------------Main--------------------------------------
int main() {
    printf("elevation angle threshold(unit:degree)(10 degree is recommended): %d\n", lim);
    printf("the order of spheric harmonic function (4 order is recommended): %d\n", order);
    printf("MDCB(multi-stations) starts running!\n");

    // Step one: read rinex files
    read_rinex(r_ipath, r_opath, &sitesInfo, &obs);

    // Step two--------------------Read SP3 files--------------------------------
    read_sp3(s_ipath, s_opath);

    // Step three------------------Read ionex files------------------------------
    read_ionex(i_ipath, i_opath, &sitesInfo, sDCB_REF);

    // Step four-------------------Ionosphere Observations-----------------------
    get_smoothed_P4(sitesInfo, (lim*pi/180), 0);


    return 0;
}

void read_rinex(const char* r_ipath, const char* r_opath, SitesInfo *psitesInfo, Obs *pobs) {
    DIR *dir;
    FILE *file;
    struct dirent *entry;
    //int counter;

    // 打开文件夹
    dir = opendir(r_ipath);
    if (dir == NULL) {
        perror("无法打开文件夹");
        exit(EXIT_FAILURE);
    }

    printf("Starting Reading Rinex Files!\n");

    // 计算文件数
    r_file_num = countFilesInDirectory(r_ipath);

    // 为 SitesInfo 结构体内的字符串数组 name 分配内存
    psitesInfo->name = (char **) malloc(r_file_num * sizeof(char *));
    for (int i = 0; i < r_file_num; i++) {
        psitesInfo->name[i] = (char *) malloc(256 * sizeof(char)); // 假设文件名最大长度为256
    }
    psitesInfo->doy = (int *) malloc(r_file_num * sizeof(int));
    for (int i = 0; i < r_file_num; i++){
        psitesInfo->doy[i] = 0;
    }
    psitesInfo->coor = (double **) malloc(r_file_num * sizeof(double *));
    for (int i = 0; i < r_file_num; i++) {
        psitesInfo->coor[i] = (double *) malloc(3 * sizeof(double)); // 假设文件名最大长度为256
    }
    for (int i = 0; i < r_file_num; i++){
        for (int j = 0; j < 3; j++){
            psitesInfo->coor[i][j] = 0;
        }
    }
    psitesInfo->RDCB_REF = (double *) malloc(r_file_num * sizeof(double));
    for (int i = 0; i < r_file_num; i++){
        psitesInfo->RDCB_REF[i] = 0;
    }

    //为obs结构体内的数组分配内存，每个变量均需要一个(2880,32)尺寸的二维double数组
    // 分配内存空间给 P1
    pobs->P1 = (double **)malloc(2880 * sizeof(double *));
    for (int i = 0; i < 2880; i++) {
        pobs->P1[i] = (double *)malloc(32 * sizeof(double));
        for (int j = 0; j < 32; j++) {
            pobs->P1[i][j] = 0.0; // 初始化为零
        }
    }
    // 分配内存空间给 P2
    pobs->P2 = (double **)malloc(2880 * sizeof(double *));
    for (int i = 0; i < 2880; i++) {
        pobs->P2[i] = (double *)malloc(32 * sizeof(double));
        for (int j = 0; j < 32; j++) {
            pobs->P2[i][j] = 0.0; // 初始化为零
        }
    }

    // 分配内存空间给 L1
    pobs->L1 = (double **)malloc(2880 * sizeof(double *));
    for (int i = 0; i < 2880; i++) {
        pobs->L1[i] = (double *)malloc(32 * sizeof(double));
        for (int j = 0; j < 32; j++) {
            pobs->L1[i][j] = 0.0; // 初始化为零
        }
    }

    // 分配内存空间给 L2
    pobs->L2 = (double **)malloc(2880 * sizeof(double *));
    for (int i = 0; i < 2880; i++) {
        pobs->L2[i] = (double *)malloc(32 * sizeof(double));
        for (int j = 0; j < 32; j++) {
            pobs->L2[i][j] = 0.0; // 初始化为零
        }
    }

    int index = 0;
    int i = 0;
    char doy[4];
    char sub_str[3];
    while ((entry = readdir(dir)) != NULL) {
        if (entry->d_name[index] != '.') {
            strcpy(psitesInfo->name[i], entry->d_name);
            strncpy(sub_str, &entry->d_name[9], 2);
            strncpy(doy, &entry->d_name[4], 3);
            sub_str[2] = '\0';
            doy[3] = '\0';  //这一步必须加入字符串标志，否则会出问题
            //printf("sub_str:%d\n", atoi(sub_str));
            //printf("doy:%d\n", atoi(doy));
            psitesInfo->doy[i] = 1000 * atoi(sub_str) + atoi(doy);
            i++;
        }
    }

    char line[256]; // 假设一行最多包含256个字符
    char filename[256];
    int obst_n; char **obst;
    double h; double m; double s; int ep = 0;
    int *loc;

    for (int i = 0; i < r_file_num; i++) {
        printf("Reading File No.%d :%s\n", i + 1, psitesInfo->name[i]);
        // 打开文件
        strcpy(filename, r_ipath);
        strcat(filename, "\\");
        strcat(filename, psitesInfo->name[i]);
        file = fopen(filename, "r");
        if (file == NULL) {
            perror("Can't Open File");
        }

        //给pobs结构体内的数组赋初值为0
        for (int j = 0; j < 2880; j++) {
            for (int k = 0; k < 32; k++) {
                pobs->P1[j][k] = 0.0;
                pobs->P2[j][k] = 0.0;
                pobs->L1[j][k] = 0.0;
                pobs->L2[j][k] = 0.0;
            }
        }

        while (fgets(line, sizeof(line), file)) {
            // -----get receivers coordinates-----------
            if (strstr(line, "APPROX POSITION XYZ") != NULL) {
                //printf("%s", line);
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
                //printf("%s", line);

                if(sscanf(line, "%d", &obst_n) == 1){
                    if (obst_n > 9){
                        printf("obst_n > 9!!!!!");
                    }
                    else{
                        obst = (char **) malloc(obst_n * sizeof(char *));
                        for (int j = 0; j < obst_n; j++) {
                            obst[j] = (char *) malloc(4 * sizeof(char)); // 假设每个观测量最大长度为4
                        }
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

                        loc = (int *) malloc(obst_n * sizeof(int));
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

//                        if (ep == 17){
//                            printf("ep: %d\n", ep);
//                        }
                    }else{
                        continue;
                    }
                    //---get satellite number----

                    //debug专用，如果line的值等于下面的字符串，则打印这行内容
                    int nsat; int *sv_G;
                    //取出line的第31到32个字符，转换成int类型，再赋值给nsat
                    if (sscanf(line + 30, "%2d", &nsat) == 1){
                    }
                    //给sv_G申请nsat个int类型的内存，并将数组中每个值赋值为0
                    sv_G = (int *) malloc(nsat * sizeof(int));
                    for (int j = 0; j < nsat; j++){
                        sv_G[j] = 0;
                    }
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
                        obs_temp = (double *) malloc(obst_n * sizeof(double));
                        //给obs_temp都赋初值为0
                        for (int k = 0; k < obst_n; k++){
                            obs_temp[k] = 0;
                        }

                        if (obst_n>5){
                            if (sv_G[j] == 0){
                                fgets(line, sizeof(line), file);
                                continue;
                            }
                            for (int k = 0; k < 5; k++){//思朗科技
                                //判断line的长度，如果大于16*(k+1)-3，则执行下面的语句
                                if (strlen(line) > 16*(k+1)-3){
                                    //strcmp(line(16*j-15:16*j-2),'              '),如果line的第15*(k+1)-15到第15*(k+1)-2个字符是空格，则执行下面的语句
                                    if (strncmp(line + 16*k, "              ", 14) == 0){
                                        continue;
                                    }
                                    //取出line的第15*(j+1)-15至第15*(j+1)-2个字符，转换成double类型，再赋值给obs_temp[j]
                                    if (sscanf(line + 16*k, "%14lf", &obs_temp[k]) == 1){
//                                        //如果obs_temp[k]的前8位数等于25789317，则打印这行内容
//                                        if (obs_temp[k] == 25789317.7094){
//                                            printf("ep: %d\n", ep);
//                                        }
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
                    }
                }
            }
        }
        // 关闭文件
        fclose(file);

        char sav_filename[256]; // 假设文件名长度不超过256个字符
        char sav_dat_filename[256]; // 假设文件名长度不超过256个字符

        //生成文件名，路径为当前路径下的RINEX_output_files文件夹
        strcpy(sav_filename, r_opath);
        strcat(sav_filename, "\\");
        strcat(sav_filename, psitesInfo->name[i]);
        strcpy(sav_dat_filename, sav_filename);
        strcat(sav_filename, ".csv");
        strcat(sav_dat_filename, ".dat");

        // Open file for writing
//        FILE *outfile = fopen(sav_filename, "w");
//        if (outfile == NULL) {
//            fprintf(stderr, "Error opening file for writing\n");
//        }

        // Write each array to file with labels
//        writeArrayToFileWithLabels(outfile, pobs->P1, "P1");
//        writeArrayToFileWithLabels(outfile, pobs->P2, "P2");
//        writeArrayToFileWithLabels(outfile, pobs->L1, "L1");
//        writeArrayToFileWithLabels(outfile, pobs->L2, "L2");

        //writeToFile(sav_dat_filename, pobs, 1);
        writeObsToFile(sav_dat_filename, pobs, 2880, 2880, 2880, 2880, 32);

        // Close file
//        fclose(outfile);

        printf("Data written to file successfully.\n");

    }
    //将obs数据以dat的格式存储
    //writeToFile("D:\\projects\\M_DCB_C\\RINEX_output_files\\sites_info.dat", &sitesInfo, r_file_num);

    printf("Step one: completings !\n");
}

void read_sp3(const char* s_ipath, const char* s_opath){
    DIR *dir;
    FILE *file;
    struct dirent *entry;
    int sp3_file_num;

    printf("Reading SP3 Files!\n");

    // 计算sp3_files文件夹中文件数
    sp3_file_num = countFilesInDirectory(s_ipath);
    //如果文件数小于3，则报错"Need at least 3 SP3 files!"并退出程序
    if (sp3_file_num < 3){
        printf("Need at least 3 SP3 files!\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < (sp3_file_num - 2); i++) {
        int index = 0;//读取文件夹中的文件的索引
        char pre_file_name[256];
        char cur_file_name[256];
        char next_file_name[256];
        int G_Week = 0;
        int Day_of_Week = 0;

        // 打开文件夹
        dir = opendir(s_ipath);
        if (dir == NULL) {
            perror("无法打开文件夹");
            exit(EXIT_FAILURE);
        }

        //给pre_file_name,cur_file_name,next_file_name都赋值
        while ((entry = readdir(dir)) != NULL) {  // 遍历文件夹中的每一个文件
            if (index == i+2) {  //前两个文件分别是"."和"..",所以从第三个文件开始读
                char filepath[256];  // 定义一个足够大的缓冲区来存储完整路径
                strcpy(filepath, s_ipath);  // 将文件夹路径拷贝到缓冲区
                strcat(filepath, "\\");  // 将文件名拷贝到缓冲区
                strcat(filepath, entry->d_name);  // 将文件名拷贝到缓冲区
                //把filepath的值赋给pre_file_name
                strcpy(pre_file_name, filepath);
            }
            if (index == i+3){
                char filepath[256];  // 定义一个足够大的缓冲区来存储完整路径
                strcpy(filepath, s_ipath);  // 将文件夹路径拷贝到缓冲区
                strcat(filepath, "\\");  // 将文件名拷贝到缓冲区
                strcat(filepath, entry->d_name);  // 将文件名拷贝到缓冲区
                strcpy(cur_file_name, filepath);

                sscanf(entry->d_name + 3, "%4d", &G_Week);
                sscanf(entry->d_name + 7, "%1d", &Day_of_Week);
            }
            if (index == i+4){
                char filepath[256];  // 定义一个足够大的缓冲区来存储完整路径
                strcpy(filepath, s_ipath);  // 将文件夹路径拷贝到缓冲区
                strcat(filepath, "\\");  // 将文件名拷贝到缓冲区
                strcat(filepath, entry->d_name);  // 将文件名拷贝到缓冲区
                strcpy(next_file_name, filepath);
            }
                index++;  // 如果不是第 i 个文件，继续查找下一个文件
        }
        closedir(dir);  // 关闭文件夹

        //读取pre_file_name sp3文件的x,y,z坐标
        //定义一个一维double数组pre_xyz，长度为3，并赋初值为0
        double ***pre_xyz; double ***cur_xyz; double ***next_xyz; double ***sate_xyz;
        pre_xyz = (double ***)malloc(3 * sizeof(double **));
        cur_xyz = (double ***)malloc(3 * sizeof(double **));
        next_xyz = (double ***)malloc(3 * sizeof(double **));
        sate_xyz = (double ***)malloc(3 * sizeof(double **));
        for (int i = 0; i < 3; i++) {
            pre_xyz[i] = (double **)malloc(96 * sizeof(double *));
            for (int j = 0; j < 96; j++) {
                pre_xyz[i][j] = (double *)malloc(32 * sizeof(double));
                for (int k = 0; k < 32; k++) {
                    pre_xyz[i][j][k] = 0.0; // 初始化为零
                }
            }
        }
        for (int i = 0; i < 3; i++) {
            cur_xyz[i] = (double **)malloc(96 * sizeof(double *));
            for (int j = 0; j < 96; j++) {
                cur_xyz[i][j] = (double *)malloc(32 * sizeof(double));
                for (int k = 0; k < 32; k++) {
                    cur_xyz[i][j][k] = 0.0; // 初始化为零
                }
            }
        }
        for (int i = 0; i < 3; i++) {
            next_xyz[i] = (double **)malloc(96 * sizeof(double *));
            for (int j = 0; j < 96; j++) {
                next_xyz[i][j] = (double *)malloc(32 * sizeof(double));
                for (int k = 0; k < 32; k++) {
                    next_xyz[i][j][k] = 0.0; // 初始化为零
                }
            }
        }
        for (int i = 0; i < 3; i++) {
            sate_xyz[i] = (double **)malloc(2880 * sizeof(double *));
            for (int j = 0; j < 2880; j++) {
                sate_xyz[i][j] = (double *)malloc(32 * sizeof(double));
                for (int k = 0; k < 32; k++) {
                    sate_xyz[i][j][k] = 0.0; // 初始化为零
                }
            }
        }

        //pre_xyz[][][],第一个[]代表x,y,z坐标，第二个[]代表行（时间），第三个[]代表列（卫星编号）
        parse_sp3(pre_file_name, pre_xyz);
        parse_sp3(cur_file_name, cur_xyz);
        parse_sp3(next_file_name, next_xyz);

        interplotation(pre_xyz, cur_xyz, next_xyz, sate_xyz);

        int* DOY = GWeek_2_DOY(G_Week, Day_of_Week);

        char sav_filename[256]; // 假设文件名长度不超过256个字符
        //生成文件名，路径为当前路径下的SP3_output_files文件夹
        strcpy(sav_filename, s_opath);
        strcat(sav_filename, "\\");
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
        save_sp3_To_Bin_File(sate_xyz, 3, 2880, 32, sav_filename);  //以二进制格式存储sp3文件
    }
}

//读取ionex文件，相关结果存放在sitesInfo.RDCB_REF和sDCB_REF结构体中
void read_ionex(const char* i_ipath, const char* i_opath, SitesInfo *psitesInfo, SDCB_REF sdcb_ref){
    int n = r_file_num;
    int new_size = 0;
    int* doys; //去重的DOY存放在doys数组中，大小为new_size
    //给doys申请n个int类型的内存，并赋初值为0
    doys = (int *) malloc(n * sizeof(int));
    for (int i = 0; i < n; i++){
        doys[i] = 0;
    }
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
    sdcb_ref.value = (double **)malloc((new_size) * sizeof(double *));
    for (int i = 0; i < new_size; i++) {
        sdcb_ref.value[i] = (double *)malloc(32 * sizeof(double));
        for (int j = 0; j < 32; j++) {
            sdcb_ref.value[i][j] = 0.0; // 初始化为零
        }
    }
    //为sdcb_ref结构体内的数组分配内存,doy需要(new_size+1)的尺寸,并赋初值为0
    sdcb_ref.doy = (int *)malloc((new_size) * sizeof(int));
    for (int i = 0; i < new_size; i++){
        sdcb_ref.doy[i] = 0;
    }

    char** ionex_filenames; int ionex_file_num = 0;

    ionex_file_num = countFilesInDirectory(i_ipath);

    //给ionex_filenames申请ionex_file_num个char*类型的内存,每个char不超过32个字符
    ionex_filenames = (char **) malloc((ionex_file_num) * sizeof(char *));
    for (int i = 0; i < ionex_file_num; i++) {
        ionex_filenames[i] = (char *) malloc(32 * sizeof(char)); // 假设每个观测量最大长度为32
    }

    get_filenames(i_ipath, ionex_filenames,ionex_file_num);

    //为DCB_rec分配内存, 大小为ionex_file_num
    DCB_rec = (double *)malloc(ionex_file_num * sizeof(double));

    //开始读文件
    for (int i= 0; i < new_size; i++){
        int index = 0; int *index2; int doy_num;
        char** sites_name;

        index2 = (int *) malloc(r_file_num * sizeof(int));
        for (int i = 0; i < r_file_num; i++){
            index2[i] = 0;
        }

        index = find_ionex_index(doys[i]);
        doy_num = find_rinex_indexs(index2, doys[i]); //在rinex文件列表中找到所有和doy[i]相等的索引号

        //给sites_name申请doy_num个char*类型的内存,每个char不超过4个字符
        sites_name = (char **) malloc(doy_num * sizeof(char *));
        for (int j = 0; j < doy_num; j++) {
            sites_name[j] = (char *) malloc(4 * sizeof(char)); // 假设每个观测量最大长度为4
        }

        //取出sites_name的值
        for (int j = 0; j < doy_num; j++){
            strncpy(sites_name[j], &psitesInfo->name[index2[j]][0],4);
            sites_name[j][4] = '\0';
        }

        //开始读ionex文件
        char ionex_file_path[256];
        strcpy(ionex_file_path, i_ipath);
        strcat(ionex_file_path, "\\");
        strcat(ionex_file_path, ionex_filenames[index]);

        read_single_ionex(ionex_file_path, i, sites_name, doy_num, sdcb_ref, DCB_rec);

        for (int j = 0; j < ionex_file_num; j++){
            sitesInfo.RDCB_REF[index2[j]] = DCB_rec[j];
        }

        sdcb_ref.doy[i] = doys[i];
    }
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
        counter ++;
//        printf("No.%d :%s\n", counter,entry->d_name); //Print file names
    }
    counter = counter - 2;

    // 关闭文件夹
    closedir(dir);

    return counter;
}

// Function to write array to file with labels
void writeArrayToFileWithLabels(FILE *file, double **array, const char* dataType) {
    // 写入列标签
    fprintf(file, "%s_Epoch", dataType); // 假设第一列为时间戳
    for (int j = 1; j <= 32; j++) {
        fprintf(file, ",%s_sv_%d", dataType, j); // 生成并写入标签
    }
    fprintf(file, "\n"); // 结束标签行

    // 写入数据
    for (int i = 0; i < 2880; i++) {
        fprintf(file, "%d", i+1); // 假设第一列为时间戳，此处简化处理
        for (int j = 0; j < 32; j++) {
            fprintf(file, ",%.4f", array[i][j]);
        }
        fprintf(file, "\n");
    }
}

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

void interplotation(double ***pre_xyz, double ***cur_xyz, double ***next_xyz, double ***sate_xyz){
    double ***cur_xyz_etd;
    cur_xyz_etd = (double ***)malloc(3 * sizeof(double **));
    for (int i = 0; i < 3; i++) {
        cur_xyz_etd[i] = (double **)malloc(105 * sizeof(double *));
        for (int j = 0; j < 105; j++) {
            cur_xyz_etd[i][j] = (double *)malloc(32 * sizeof(double));
            for (int k = 0; k < 32; k++) {
                cur_xyz_etd[i][j][k] = 0.0; // 初始化为零
            }
        }
    }

    extend_matrix(0, pre_xyz, cur_xyz, next_xyz, cur_xyz_etd);//x_extend
    extend_matrix(1, pre_xyz, cur_xyz, next_xyz, cur_xyz_etd);//y_extend
    extend_matrix(2, pre_xyz, cur_xyz, next_xyz, cur_xyz_etd);//z_extend

    //创建一个整数数组m_t，起始值为-120，结束值为3000，数组长度为105
    int m_t[105];
    for (int i = 0; i < 105; i++){
        m_t[i] = -120 + 30*i;
    }

    for (int i = 0; i < 32; i++){
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
        }
    }

}

void extend_matrix(int dimension, double ***pre_xyz, double ***cur_xyz, double ***next_xyz, double ***xyz_etd){
    //从矩阵pre_xyz[0]中选取第93到第96行的所有列，形成一个子矩阵。将子矩阵与矩阵cur_xyz[0]进行垂直拼接，即将子矩阵放在cur_xyz[0]的上方。将矩阵next_xyz[0]中的第1到第5行的所有列与前一步得到的结果进行垂直拼接，即将next_xyz[0]的前5行放在前一步得到的结果的下方。
    // 先填充pre_xyz的部分
    for (int j = 0; j < 4; j++){
        for (int k = 0; k < 32; k++) {
            xyz_etd[dimension][j][k] = pre_xyz[dimension][j+92][k];
        }
    }

    // 填充cur_xyz的部分
    for (int j = 4; j < 100; j++) {
        for (int k = 0; k < 32; k++) {
            xyz_etd[dimension][j][k] = cur_xyz[dimension][j-4][k];
        }
    }

    // 填充next_xyz的部分
    for (int j = 100; j < 105; j++) {
        for (int k = 0; k < 32; k++) {
            xyz_etd[dimension][j][k] = next_xyz[dimension][j-100][k];
        }
    }
}

// 声明Lagrange插值函数
double* interp_lag(double* x, double* y, double* x0) {
    int n = 10; // 假设x中有10个点
    int n0 = 30;
    double* y0 = (double*)malloc(n0 * sizeof(double)); // 动态分配返回数组
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

int* GWeek_2_DOY(int G_Week, int Day_of_Week){
    int *DOY = (int *) malloc(2 * sizeof(int));
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
    for (int svn = 1; svn <= 32; svn++) {
        fprintf(file, "%d", svn);
        if (svn < 32) fprintf(file, ",");
        else fprintf(file, "\n");
    }

    // 分别按照x、y、z的顺序写入数据
    for (int dim = 0; dim < 3; dim++) { // 0: x, 1: y, 2: z
        for (int epoch = 0; epoch < 2880; epoch++) {
            fprintf(file, "Epoch %d,", epoch + 1);
            for (int svn = 0; svn < 32; svn++) {
                fprintf(file, "%lf", sate_xyz[dim][epoch][svn]);
                if (svn < 31) fprintf(file, ",");
                else fprintf(file, "\n");
            }
        }
        // 在不同维度数据间添加空行，增加可读性，可根据个人喜好决定是否保留
        fprintf(file, "\n");
    }

    fclose(file); // 关闭文件
    printf("SP3 data successfully saved to %s\n", filename);
}

int find_ionex_index(int doy){
    int index = 0;

    //读取i_ipath文件夹中的文件
    DIR *dir;
    FILE *file;
    struct dirent *entry;

    // 打开文件夹
    dir = opendir(i_ipath);
    if (dir == NULL) {
        perror("无法打开文件夹");
        exit(EXIT_FAILURE);
    }

    while ((entry = readdir(dir)) != NULL) {
        int number1 = 0; int number2 = 0; int doy_temp = 0;
        //如果entry->d_name的第一个字符是"i"，则跳入下一个循环
        if (entry->d_name[0] == '.' || entry->d_name == ".."){
            continue;
        }
        //第5到第7位的int类型数字1
        sscanf(entry->d_name + 4, "%3d", &number1);
        //第10到第11位的int类型数字2
        sscanf(entry->d_name + 9, "%2d", &number2);
        doy_temp = number1 + number2*1000;
        if (doy_temp == doy){
            return index;
            break;
        }
        index++;
    }

    perror("Can't find the ionex file");
    return -1;
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

void get_filenames(const char* folderPath, char** filenames, int file_num) {
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
        if (entry->d_name[0] == '.' || entry->d_name == ".."){
            continue;
        }
        strcpy(filenames[counter], entry->d_name);
        counter++;
    }

    // 关闭文件夹
    closedir(dir);
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
    //依次读取每天的rinex参数和sp3参数
    for (int i=0; i < r_file_num; i++){
        //读取RINEX_output_files文件夹中的第一个文件
        char rinex_file_path[256];
        Obs obs_temp;
        strcpy(rinex_file_path, r_opath);
        strcat(rinex_file_path, "\\");
        strcat(rinex_file_path, sitesInfo.name[i]);
        strcat(rinex_file_path, ".dat");
        readObsFromFile(rinex_file_path, &obs_temp);

        char* site; int doy;
        site = (char *) malloc(64 * sizeof(char));
        sscanf(sitesInfo.name[i], "%4s", site);
        doy = sitesInfo.doy[i];
        int index = 0;

        //找到去重后的rinex文件列表中与site对应的index
        int new_size = 0;
        char** sites; //去重的DOY存放在doys数组中，大小为new_size
        //给doys申请n个int类型的内存，并赋初值为0
        sites = (char **) malloc(r_file_num * sizeof(char));
        for (int j = 0; j < r_file_num; j++){
            sites[j] = (char *) malloc(4 * sizeof(char)); // 假设每个观测量最大长度为4
        }
        //遍历r_file的doy
        sscanf(sitesInfo.name[0], "%4s", sites[0]);
        new_size++;
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
        }

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
        strcat(sp3_load_file_path, "\\");
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
        sate_xyz = (double ***)malloc(3 * sizeof(double **));
        for (int j = 0; j < 3; j++) {
            sate_xyz[j] = (double **)malloc(2880 * sizeof(double *));
            for (int k = 0; k < 2880; k++) {
                sate_xyz[j][k] = (double *)malloc(32 * sizeof(double));
                for (int l = 0; l < 32; l++) {
                    sate_xyz[j][k][l] = 0.0; // 初始化为零
                }
            }
        }

        load_Sp3_From_Bin_File(&sate_xyz, 3, 2880, 32, sp3_load_file_path);

        cutobs(&obs_temp, sate_xyz, sx, sy, sz, z_threshold);//筛选后的obs值存在obs_temp中



    }
}

// 将结构体数组写入文件
void writeToFile(char* filename, Obs* obs, int num_sites) {
    FILE* fp = fopen(filename, "wb");
    if (fp == NULL) {
        printf("无法打开文件 %s\n", filename);
        return;
    }

    // 写入结构体数组数据到文件中
    fwrite(obs, sizeof(obs), num_sites, fp);

    fclose(fp);
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
    obs->P1 = (double**)malloc(p1Rows * sizeof(double*));
    for (int i = 0; i < p1Rows; i++) {
        obs->P1[i] = (double*)malloc(cols * sizeof(double));
        fread(obs->P1[i], sizeof(double), cols, file);
    }
    // 重复上述过程读取P2, L1, L2
    obs->P2 = (double**)malloc(p2Rows * sizeof(double*));
    for (int i = 0; i < p2Rows; i++) {
        obs->P2[i] = (double*)malloc(cols * sizeof(double));
        fread(obs->P2[i], sizeof(double), cols, file);
    }
    obs->L1 = (double**)malloc(l1Rows * sizeof(double*));
    for (int i = 0; i < l1Rows; i++) {
        obs->L1[i] = (double*)malloc(cols * sizeof(double));
        fread(obs->L1[i], sizeof(double), cols, file);
    }
    obs->L2 = (double**)malloc(l2Rows * sizeof(double*));
    for (int i = 0; i < l2Rows; i++) {
        obs->L2[i] = (double*)malloc(cols * sizeof(double));
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

void cutobs(Obs *obs_temp, double*** sate_xyz, double sx, double sy, double sz, int lim){
    double** x = sate_xyz[0];
    double** y = sate_xyz[1];
    double** z = sate_xyz[2];

    for (int i = 0; i < 32; i++){
        for (int j = 0; j < 2880; j++){
            if (obs_temp->L1[j][i] == 0 || obs_temp->L2[j][i] == 0 || obs_temp->P1[j][i] == 0 || obs_temp->P2[j][i] == 0){
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
    double* EA = (double *) malloc(2 * sizeof(double));
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
    double* BLH = (double *) malloc(3 * sizeof(double));
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
    double** L6 = (double **)malloc(2880 * sizeof(double *));
    for (int i = 0; i < 2880; i++) {
        L6[i] = (double *)malloc(32 * sizeof(double));
        for (int j = 0; j < 32; j++) {
            L6[i][j] = 0.0; // 初始化为零
        }
    }

    double** Li = (double **)malloc(2880 * sizeof(double *));
    for (int i = 0; i < 2880; i++) {
        Li[i] = (double *)malloc(32 * sizeof(double));
        for (int j = 0; j < 32; j++) {
            Li[i][j] = 0.0; // 初始化为零
        }
    }

    double** Nw = (double **)malloc(2880 * sizeof(double *));
    for (int i = 0; i < 2880; i++) {
        Nw[i] = (double *)malloc(32 * sizeof(double));
        for (int j = 0; j < 32; j++) {
            Nw[i][j] = 0.0; // 初始化为零
        }
    }

    double** P4 = (double **)malloc(2880 * sizeof(double *));
    for (int i = 0; i < 2880; i++) {
        P4[i] = (double *)malloc(32 * sizeof(double));
        for (int j = 0; j < 32; j++) {
            P4[i][j] = 0.0; // 初始化为零
        }
    }

    double** L4 = (double **)malloc(2880 * sizeof(double *));
    for (int i = 0; i < 2880; i++) {
        L4[i] = (double *)malloc(32 * sizeof(double));
        for (int j = 0; j < 32; j++) {
            L4[i][j] = 0.0; // 初始化为零
        }
    }

    // Wide Lane Wavelength
    double lambda_w = c/(f1-f2);
    // MW observable
    //L6 = lambda_w*(obs_temp->L1 - obs_temp->L2) - (f1*obs.P1+f2*obs.P2)/(f1+f2);
    for (int i = 0; i < 2880; i++){
        for (int j = 0; j < 32; j++){
            L6[i][j] = lambda_w*(obs_temp->L1[i][j] - obs_temp->L2[i][j]) - (f1*obs_temp->P1[i][j]+f2*obs_temp->P2[i][j])/(f1+f2);
        }
    }
    //Li = obs.L1 - f1*obs.L2/f2;
    for (int i = 0; i < 2880; i++){
        for (int j = 0; j < 32; j++){
            Li[i][j] = obs_temp->L1[i][j] - f1*obs_temp->L2[i][j]/f2;
        }
    }

    //Nw = L6/lambda_w;
    for (int i = 0; i < 2880; i++){
        for (int j = 0; j < 32; j++){
            Nw[i][j] = L6[i][j]/lambda_w;
        }
    }

    for (int prn = 0; prn < 32; prn++){
        int** arc; int arc_n = 0; int aaa = 0;

        //------divide arc---------------------------
        arc = Get_arc(L6, prn);//arc是一个n行2列的数组，每一行代表一个arc的起始和终止的index
        //arc_n是arc的行数
        arc_n = sizeof(arc) / sizeof(arc[0]);
        //aaa是arc的列数
        aaa = sizeof(arc[0]) / sizeof(arc[0][0]);

        //----delete arc less than 10 epoches-------
        int* arc_d; int index_of_arc_d = 0;
        arc_d = (int *) malloc(arc_n * sizeof(int));//用于记录epoch<10在arc_n中的索引

        for (int j = 0; j < arc_n; j++){
            int n_epoch = arc[j][1] - arc[j][0];
            if (n_epoch<10){
                for (int k = arc[j][0]; k < arc[j][1]; k++){
                    obs_temp->P1[k][prn] = 0; obs_temp->P2[k][prn] = 0; obs_temp->L1[k][prn] = 0; obs_temp->L2[k][prn] = 0;
                    L6[k][prn] = 0; Li[k][prn] = 0; Nw[k][prn] = 0;
                }
                arc_d[index_of_arc_d] = j;
                index_of_arc_d++;
            }
        }

        int d_size = sizeof(arc_d) / sizeof(arc_d[0]);//需要删除的行数量
        int new_size;

        int** arc_final = removeRows(arc, arc_n, arc_d, d_size, &new_size);//删除对应行后的arc数组，之后都处理arc_final数组

        //----mw detect cycle slip------------------
        int arc_final_n = sizeof(arc_final) / sizeof(arc_final[0]);
        int arc_final_aaa = sizeof(arc_final[0]) / sizeof(arc_final[0][0]);
        int j = 0;

        while (j<arc_final_n){
            //----first epoch check----------
            int e = arc_final[j][0];
            while (1){
                if (e+1 == arc_final[j][1] || e == arc_final[j][1]){
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
                    arc_final[j][0] = e;
                }else{
                    arc_final[j][0] = e;
                    break;
                }
            }

            //----detect------------------
            if (arc_final[j][1] - arc_final[j][0]<10){
                for (int k = arc_final[j][0]; k < arc_final[j][1]; k++){
                    obs_temp->P1[k][prn] = 0; obs_temp->P2[k][prn] = 0; obs_temp->L1[k][prn] = 0; obs_temp->L2[k][prn] = 0;
                    L6[k][prn] = 0; Li[k][prn] = 0; Nw[k][prn] = 0;
                }
                //删除掉arc_final[j]这一行
                deleteRowAndReplace(&arc_final, &arc_final_n, j);
                arc_final_n--;
                continue;
            }

            double* ave_N = (double *) malloc((arc_final[j][1] - arc_final[j][0] + 1) * sizeof(double));
            double* sigma = (double *) malloc((arc_final[j][1] - arc_final[j][0] + 1) * sizeof(double));
            double* sigma2 = (double *) malloc((arc_final[j][1] - arc_final[j][0] + 1) * sizeof(double));
            ave_N[0] = Nw[arc_final[j][0]][prn];
            sigma[0] = 0; sigma2[0] = 0;
            int count = 2;

            //----------------------check epoch k+1


        }

        //--------smoothing-------------------------


        //--------remove bad P4---------------------


    }
}

int** Get_arc(double** L6, int prn){
    int* arc; int index_of_arc = 0;
    arc = (int *) malloc(2880 * sizeof(int));
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
            arc[index_of_arc] = i;
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

    //遍历arc数组，如果除了第1个值的其他值不等于0，则计数器加1，以计算arc数组的长度
    int counter = 0;
    for (int i = 0; i < 2880; i++){
        if (i ==0 && arc[0] == 0 && arc[1] != 0){
            counter++;
        }
        if (arc[i] != 0){
            counter++;
        }
    }

    //创建一个新的数组arc_temp，长度为counter
    int* arc_temp;
    arc_temp = (int *) malloc(counter * sizeof(int));
    for (int i = 0; i < counter; i++){
        arc_temp[i] = arc[i];
    }

    int size = sizeof(arc_temp) / sizeof(arc_temp[0]);
    int rows = size / 2;

    // 动态分配二维数组的内存
    int ** result;
    //申请一个rows行2列的内存
    result = (int **)malloc(rows * sizeof(int *));
    for (int i = 0; i < rows; i++) {
        result[i] = (int *)malloc(2 * sizeof(int));
    }

    convertTo2D(arc_temp, size, result);

    return result;
}

void convertTo2D(int *arr, int size, int** result){
    // 检查元素数量是否为偶数
    if (size % 2 != 0) {
        printf("Error: The number of elements in the array must be even.\n");
        exit(1);
    }

    // 将一维数组的元素复制到二维数组中
    for (int i = 0; i < size / 2; i++) {
        result[i][0] = arr[2 * i];
        result[i][1] = arr[2 * i + 1];
    }
}

int** removeRows(int** arc, int n, int* arc_d, int d_size, int* new_size) {
    // 计算新数组的行数
    *new_size = n - d_size;

    // 分配新数组的内存
    int** arc_final = (int**)malloc(*new_size * sizeof(int*));
    for (int i = 0; i < *new_size; i++) {
        arc_final[i] = (int*)malloc(2 * sizeof(int)); // 假设二维数组的第二维大小固定为2
    }

    int arc_d_index = 0, arc_final_index = 0;
    for (int i = 0; i < n; i++) {
        if (arc_d_index < d_size && i == arc_d[arc_d_index]) {
            // 如果当前行索引在删除列表中，跳过
            arc_d_index++;
        } else {
            // 否则，复制到新数组中
            arc_final[arc_final_index][0] = arc[i][0];
            arc_final[arc_final_index][1] = arc[i][1];
            arc_final_index++;
        }
    }

    return arc_final;
}

void deleteRowAndReplace(int ***arc, int *n, int row) {
    int newN = *n - 1; // 新数组的行数
    int **tempArc = (int **)malloc(newN * sizeof(int *));
    for (int i = 0; i < newN; i++) {
        tempArc[i] = (int *)malloc(2 * sizeof(int));
    }

    // 复制除了要删除的行之外的所有行到临时数组
    for (int i = 0, j = 0; i < *n; i++) {
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