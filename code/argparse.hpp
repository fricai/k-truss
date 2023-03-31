#pragma once

#include <string>
#include <cstring>
#include <vector>
#include <iostream>

struct arguments_t {
    int taskid;
    std::string inputpath;
    std::string headerpath;
    std::string outputpath;
    bool verbose;
    int startk;
    int endk;
    int p;
};

#ifndef __clang__

#include <argp.h>

const char* argp_program_version = "Unroll v.0.1";
const char* argp_program_bug_address = "cs1200837_cs1200869";
static char doc[] = "fast k-truss decomposition on graphs";
static char args_doc[] = "";
static argp_option options[] = {
    {"taskid", 't', "taskid", 0, "Task ID (defaults to 1)"},
    {"inputpath", 'i', "inputpath", 0, "Path to graph input file"},
    {"headerpath", 'h', "headerpath", 0, "Path to header file"},
    {"outputpath", 'o', "outputpath", 0, "Path to output file"},
    {"verbose", 'v', "verbose", 0, "Verbose output (defaults to 0)"},
    {"startk", 's', "startk", 0, "Start of range"},
    {"endk", 'e', "endk", 0, "End of range"},
    {"p", 'p', "p", 0, "Influencer degree"},
    {0}};

static error_t parse_opt(int key, char* arg, struct argp_state* state) {
    arguments_t* arguments = (arguments_t*)state->input;
    switch (key) {
        case 't':
            arguments->taskid = arg[0] - '0';
            break;
        case 'i':
            arguments->inputpath = std::string(arg);
            break;
        case 'h':
            arguments->headerpath = std::string(arg);
            break;
        case 'o':
            arguments->outputpath = std::string(arg);
            break;
        case 'v':
            arguments->verbose = (arg[0] == '1');
            break;
        case 's':
            arguments->startk = std::stoi(arg);
            break;
        case 'e':
            arguments->endk = std::stoi(arg);
            break;
        case 'p':
            arguments->p = std::stoi(arg);
            break;
        case ARGP_KEY_ARG:
            return 0;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static argp argp = {options, parse_opt, args_doc, doc, 0, 0};

#else

static int argp = 0;

void argp_parse(int* argp, int argc, char** argv, int opt1, int opt2, arguments_t *args) {

    // do some argparsing here
    std::vector<bool> proc_args(argc, false);

    for (int i=0; i<argc; i++) {
        if (proc_args[i]) continue;
        if (i == 0) {
            // program name
            proc_args[0] = true;
            continue;
        }
        if (strcmp(argv[i],"--inputpath") == 0) {
            args->inputpath = std::string(argv[i+1]);
            proc_args[i] = proc_args[i+1] = true;
            continue;
        }
        if (strcmp(argv[i],"--headerpath") == 0) {
            args->headerpath = std::string(argv[i+1]);
            proc_args[i] = proc_args[i+1] = true;
            continue;
        }
        if (strcmp(argv[i],"--startk") == 0) {
            args->startk = std::stoi(argv[i+1]);
            proc_args[i] = proc_args[i+1] = true;
            continue;
        }
        if (strcmp(argv[i],"--endk") == 0) {
            args->endk = std::stoi(argv[i+1]);
            proc_args[i] = proc_args[i+1] = true;
            continue;
        }
        if (strcmp(argv[i],"--outputpath") == 0) {
            args->outputpath = std::string(argv[i+1]);
            proc_args[i] = proc_args[i+1] = true;
            continue;
        }
        if (strcmp(argv[i],"--taskid") == 0) {
            args->outputpath = std::stoi(argv[i+1]);
            proc_args[i] = proc_args[i+1] = true;
            continue;
        }
        if (strcmp(argv[i],"--verbose") == 0) {
            args->verbose = std::stoi(argv[i+1]);
            proc_args[i] = proc_args[i+1] = true;
            continue;
        }
        if (strcmp(argv[i],"--p") == 0) {
            args->p = std::stoi(argv[i+1]);
            proc_args[i] = proc_args[i+1] = true;
            continue;
        }
    }

}

#endif

arguments_t args;
