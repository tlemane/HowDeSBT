#ifndef cmd_compress_bf_H
#define cmd_compress_bf_H

#include <string>
#include <cstdlib>
#include <cstdint>
#include <iostream>

#include "commands.h"

class CompressBFCommand: public Command
	{
public:
	CompressBFCommand(const std::string& name): Command(name) {}
	virtual ~CompressBFCommand() {}
	virtual void short_description (std::ostream& s);
	virtual void usage (std::ostream& s, const std::string& message="");
	virtual void debug_help (std::ostream& s);
	virtual void parse (int _argc, char** _argv);
	virtual int execute (void);
	virtual std::string process_bloom_filter (const std::string& filename);

	std::vector<std::string> bfFilenames;
	std::string listFilename;
	std::string inTreeFilename;
	std::string outTreeFilename;
	std::uint32_t compressor;
	};

#endif // cmd_compress_bf_H