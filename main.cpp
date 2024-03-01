#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <math.h>

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
void r_ionex(const char* ionex_file, int i);

// -----------------------------Main--------------------------------------
int main() {
    printf("elevation angle threshold(unit:degree)(10 degree is recommended): %d\n", lim);
    printf("the order of spheric harmonic function (4 order is recommended): %d\n", order);
    printf("MDCB(multi-stations) starts running!\n");

    // Step one: read rinex files
    read_rinex(r_ipath, r_opath, &sitesInfo, &obs);

    // Step two--------------------Read SP3 files--------------------------------
    //read_sp3(s_ipath, s_opath);

    // Step three------------------Read ionex files------------------------------
    read_ionex(i_ipath, i_opath, &sitesInfo, sDCB_REF);

    // 释放 sitesInfo 结构体内存
    //free_usage(sitesInfo, len);

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

        //生成文件名，路径为当前路径下的RINEX_output_files文件夹
        strcpy(sav_filename, r_opath);
        strcat(sav_filename, "\\");
        strcat(sav_filename, psitesInfo->name[i]);
        strcat(sav_filename, ".csv");

        // Open file for writing
        FILE *outfile = fopen(sav_filename, "w");
        if (outfile == NULL) {
            fprintf(stderr, "Error opening file for writing\n");
        }

        // Write each array to file with labels
        writeArrayToFileWithLabels(outfile, pobs->P1, "P1");
        writeArrayToFileWithLabels(outfile, pobs->P2, "P2");
        writeArrayToFileWithLabels(outfile, pobs->L1, "L1");
        writeArrayToFileWithLabels(outfile, pobs->L2, "L2");

        // Close file
        fclose(outfile);

        printf("Data written to file successfully.\n");

    }
    printf("Step one: completing !\n");
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
        strcat(sav_filename, ".csv");

        saveSp3ToCSV(sate_xyz, sav_filename);
    }
}

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

    //为sdcb_ref结构体内的数组分配内存,value需要((new_size+1)*32)的尺寸
    sdcb_ref.value = (double **)malloc((new_size+1) * sizeof(double *));
    for (int i = 0; i < new_size+1; i++) {
        sdcb_ref.value[i] = (double *)malloc(32 * sizeof(double));
        for (int j = 0; j < 32; j++) {
            sdcb_ref.value[i][j] = 0.0; // 初始化为零
        }
    }
    //为sdcb_ref结构体内的数组分配内存,doy需要(new_size+1)的尺寸,并赋初值为0
    sdcb_ref.doy = (int *)malloc((new_size+1) * sizeof(int));
    for (int i = 0; i < new_size+1; i++){
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
    for (int i= 0; i < new_size+1; i++){
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

        r_ionex(ionex_file_path, i);

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
            if (sscanf(line+18, "%14lf", &xyz[2][ep-1][sv-1]) == 1){
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

void r_ionex(const char* ionex_file, int i) {
    int flag = 0;

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
                if (strlen(line) > 75 && (strncmp(line+3,"G",1)||strncmp(line+3," ",1)) && strncmp(line + 60, "PRN / BIAS / RMS", 16)){
                    sscanf(line + 4, "%2d", &prn);
                    sscanf(line + 9, "%7lf", &sDCB_REF.value[i]);
                    continue;
                }
                //--receivers' DCB
                if (strlen(line) > 10 && (strncmp(line+3,"G",1)||strncmp(line+3," ",1)) && strncmp(line + 60, "STATION / BIAS / RMS", 20)){
                    
                }
            }

        }
    }
}