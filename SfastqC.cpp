#include <cstdlib>
#include <unistd.h>

#include "SfastqCManager.h"
#include "pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"

using namespace std;
using namespace PgTools;

int main(int argc, char *argv[])
{
    int opt; // current option
    PgRCManager* pgRC = new PgRCManager();
    bool expectedPairFile = false;
    bool srcFilePresent = false;
    bool pairFilePresent = false;
    bool compressionParamPresent = false;
    bool decompressMode = false;

#ifndef DEVELOPER_BUILD
    NullBuffer null_buffer;
    std::ostream null_stream(&null_buffer);
    logout = &null_stream;
#endif

#ifdef DEVELOPER_BUILD
    while ((opt = getopt(argc, argv, "c:i:q:g:s:M:p:l:B:E:doSIrNVvtaA?")) != -1) {
        char* valPtr;
#else
    while ((opt = getopt(argc, argv, "c:i:q:g:s:M:p:do?")) != -1) {
#endif
        switch (opt) {
            case 'c':
                compressionParamPresent = true;
                pgRC->setCompressionLevel(atoi(optarg));
                break;
            case 'i':
                pgRC->setSrcFastqFile(optarg);
                srcFilePresent = true;
                if (optind < (argc - 1) && argv[optind][0] != '-') {
                    pgRC->setPairFastqFile(argv[optind++]);
                    pairFilePresent = true;
                }
                break;
            case 'o':
                compressionParamPresent = true;
                pgRC->setPreserveOrderMode();
                break;
            case 'd':
                decompressMode = true;
                break;
#ifdef DEVELOPER_BUILD
                valPtr = optarg + 1;
                switch (*optarg) {
                    case 'd':case 'i':case 'c':
                        if(*valPtr == 's') {
                            valPtr++;
                            pgRC->setMatchingMode(toupper(*optarg));
                        } else
                            pgRC->setMatchingMode(*optarg);
                        break;
                    default: valPtr--;
                }
                pgRC->setReadSeedLength(atoi(valPtr));
#else
                pgRC->setReadSeedLength(atoi(optarg));
#endif
                break;
#ifdef DEVELOPER_BUILD
            case 'l':
                compressionParamPresent = true;
                valPtr = optarg + 1;
                switch (*optarg) {
                    case 'd':case 'i':case 'c':
                        if(*valPtr == 's') {
                            valPtr++;
                            pgRC->setPreMatchingMode(toupper(*optarg));
                        } else
                            pgRC->setPreMatchingMode(*optarg);
                        break;
                    default: valPtr--;
                }
                pgRC->setPreReadsExactMatchingChars(atoi(valPtr));
                break;
            case 't':
                compressionParamPresent = true;
                plainTextWriteMode = true;
                break;
#endif
            case '?':
            default: /* '?' */
                fprintf(stderr, "SfastqCManager\n\n");
                fprintf(stderr, "Usage: %s [-c compressionLevel] [-i fastqSrcFile [pairSrcFile]] [-o] [-d] "
                                "archiveName\n\n", argv[0]);
                fprintf(stderr, "-c compression levels: 1 - fast; 2 - default; 3 - max\n");
                fprintf(stderr, "-d decompression mode\n");
                fprintf(stderr, "-o preserve original read order information\n\n");
#ifdef DEVELOPER_BUILD
                                "[matchingMode]"
#endif
                                "lengthOfReadSeedPartForReadsAlignmentPhase]\n"
#ifdef DEVELOPER_BUILD
              
#endif
                fprintf(stderr, "The order of all selected options is arbitrary.\n\n");
                exit(EXIT_FAILURE);
        }
    }
    if (optind > (argc - 1) || optind < (argc - 1)) {
        fprintf(stderr, "%s: Expected 1 argument after options (found %d)\n", argv[0], argc - optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if (decompressMode && compressionParamPresent) {
        fprintf(stderr, "Cannot use compression options in decompression mode.\n");
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if (!srcFilePresent && !decompressMode) {
        fprintf(stderr, "Input file(s) not specified.\n");
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        exit(EXIT_FAILURE);
    }
#ifdef DEVELOPER_BUILD
    if (expectedPairFile && !pairFilePresent) {
        fprintf(stderr, "Cannot use -r or -I option without specifying a pair file.\n");
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        exit(EXIT_FAILURE);
    }
#endif

    pgRC->setPgRCFileName(argv[optind++]);

    if (decompressMode)
        pgRC->decompressPgRC();
    else
        pgRC->executePgRCChain();

    delete(pgRC);

    exit(EXIT_SUCCESS);
}