#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <typeinfo>

// constant info
const std::string PROGNAME = "compare_order";
const std::string VERSION = "1.0.5";
const std::string USAGE = PROGNAME + std::string(" v") + VERSION + std::string("\nUsage:\ncompare_order --file1 in1.fq --file2 in2.fq --format fastq\ncompare_order --file1 in1.txt --file2 in2.txt --format scarf\n");

// globals
std::string file1;
std::string file2;
std::string format;

void error( std::string const& err )
{
	std::cerr << std::string("Error: ") << err << std::endl;
	std::exit(1);
}

template <class T>
std::string convertToString( const T& t )
{
	std::stringstream ss;
	if ( (ss << t).fail() )
	{
		std::string s = typeid(t).name();
		error("Program died converting " + s + " to string");
	}
	return ss.str();
};

void parseOptions( int argCount, char *commands[] )
{
	if (argCount != 7)
	{
		error( "Incorrect number of parameters used!\n\n" + USAGE );
	}

	int cmdIndex = 0;

	while (++cmdIndex < argCount)
	{
		if (std::strcmp(commands[cmdIndex], "--format") == 0)
		{
			if (++cmdIndex < argCount)
			{
				if (std::strcmp(commands[cmdIndex], "fastq") == 0 || std::strcmp(commands[cmdIndex], "scarf") == 0)
				{
					format = commands[cmdIndex];
				}
				else
				{
					error("--format [" + std::string(commands[cmdIndex]) + "] not recognized!");
				}
			}
			else
			{
				error( "--format used but no format specified!" );
			}
		}
		else if (std::strcmp(commands[cmdIndex], "--file1") == 0)
		{
			if (++cmdIndex < argCount)
			{
				file1 = commands[cmdIndex];
			}
			else
			{
				error( "--file1 used with no value specified!");
			}
		}
		else if (std::strcmp(commands[cmdIndex], "--file2") == 0)
		{
			if (++cmdIndex < argCount)
			{
				file2 = commands[cmdIndex];
			}
			else
			{
				error( "--file2 used with no value specified!");
			}
		}
		else
		{
			error("Command-line value [" + std::string(commands[cmdIndex]) + "] not recognized!");
		}
	}

	if (format.empty() || file1.empty() || file2.empty())
	{
		error("Please specify --file1, --file2 and --format");
	}

	if (file1 == file2)
	{
		error("--file1 and --file2 have the same file name! Nothing to do.");
	}

	std::cout << PROGNAME << " v" << VERSION << std::endl;
	std::cout << std::string("--file1: ") << file1 << std::endl;
	std::cout << std::string("--file2: ") << file2 << std::endl;
	std::cout << std::string("--format: ") << format << std::endl;
}

void cleanHeader( std::string &s, unsigned long long int count )
{
	size_t search_pos = s.find_last_of( std::string( "#" ) );
	
	if ( search_pos == std::string::npos )
	{
		search_pos = s.find_last_of( std::string( " " ) );
		
		if ( search_pos == std::string::npos )
		{
			search_pos = s.find_last_of( std::string( "/" ) );
			
			if ( search_pos == std::string::npos )
			{
				error( "No \'#\' or \' \' character in read 1 at line " + convertToString( count ) );
			}
		}
	}

	s.erase( search_pos );
}

int main( int argc, char *argv[])
{
	clock_t tStart = std::clock();
	
	unsigned long long int lineCount = 0;
	unsigned long long int readCount = 0;
	
	parseOptions( argc, argv );
	
	std::ifstream IFH1( file1.c_str(), std::ios::in );
	std::ifstream IFH2( file2.c_str(), std::ios::in );

	if (!IFH1.good())
	{
		error( "Couldn't open file 1 " + file1 ) ;
	}

	if (!IFH2.good())
	{
		error( "Couldn't open file 2 " + file2 );
	}

	std::string buffer1;
	std::string buffer2;

	if (format == std::string("scarf"))
	{
		while (++lineCount && ++readCount)
		{
			if (IFH1.eof() && IFH2.eof())
			{
				break;
			}
			else if (IFH1.eof())
			{
				error( "File 1 ended but file 2 did not! Line " + convertToString( lineCount ) + "." );
			}
			else if (IFH2.eof())
			{
				error( "File 2 ended but file 1 did not! Line " + convertToString( lineCount ) + "." );
			}
			
			std::getline( IFH1, buffer1 );
			std::getline( IFH2, buffer2 );

			if (buffer1.empty() && buffer2.empty())
			{
				continue;
			}
			else if (buffer1.empty())
			{
				error( "File 1 did not read input at line " + convertToString( lineCount ) );
			}
			else if (buffer2.empty())
			{
				error( "File 2 did not read input at line " + convertToString( lineCount ) );
			}

			cleanHeader( buffer1, lineCount );
			cleanHeader( buffer2, lineCount );
			
			if (buffer1 != buffer2)
			{
				error( "Read 1 ID != Read 2 ID at line " + convertToString( lineCount ) + "\nFile 1: " + buffer1 + "\nFile 2: " + buffer2 + "\n" );
			}
			
			buffer1.clear();
			buffer2.clear();
		}
	}
	else if (format == std::string("fastq"))
	{
		short int fqLine = 0;
		
		while (++fqLine && ++lineCount)
		{
			if (IFH1.eof() && IFH2.eof())
			{
				break;
			}
			else if (IFH1.eof())
			{
				error( "File 1 ended but file 2 did not! Line " + convertToString( lineCount ) + "." );
			}
			else if (IFH2.eof())
			{
				error( "File 2 ended but file 1 did not! Line " + convertToString( lineCount ) + "." );
			}
			
			std::getline( IFH1, buffer1 );
			std::getline( IFH2, buffer2 );

			if (buffer1.empty() && buffer2.empty())
			{
				continue;
			}
			else if (buffer1.empty())
			{
				error( "File 1 did not read input at line " + convertToString( lineCount ) );
			}
			else if (buffer2.empty())
			{
				error( "File 2 did not read input at line " + convertToString( lineCount ) );
			}

			if (fqLine == 1)
			{
				cleanHeader( buffer1, lineCount );
				cleanHeader( buffer2, lineCount );

				if (buffer1 != buffer2)
				{
					error( "Read 1 ID != Read 2 ID at line " + convertToString( lineCount ) + "\nFile 1: " + buffer1 + "\nFile 2: " + buffer2 + "\n" );
				}
			}
			else if (fqLine == 4)
			{
				fqLine = 0;
				++readCount;
			}
			
			buffer1.clear();
			buffer2.clear();
		}
	}
	
	IFH1.close();
	IFH2.close();
	
	clock_t tEnd = std::clock();

	std::cout << std::string("All reads in sync. ") << readCount << std::string(" total reads processed in each file.") << std::endl;
	std::cout << std::string("Analysis complete in ") << std::setprecision(3) << (static_cast<float>(tEnd)-static_cast<float>(tStart)) / CLOCKS_PER_SEC << std::string(" seconds") << std::endl;
	
	return 0;
}
