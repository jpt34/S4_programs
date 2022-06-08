#pragma once

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <string>
#include <fstream>
#include <jlibrary.h>


#define ZINPUT(input, inputfile) zifstream input(inputfile)

#define ZOUTPUT(output, outputfile) zofstream output(outputfile)


namespace
{
	void CreateInputFilter(const std::string& inputfile, boost::iostreams::filtering_istream& input)
	{	
		if(inputfile.substr(inputfile.size()-3) == ".gz")
		{	
            input.push(boost::iostreams::gzip_decompressor());
		}

        if(inputfile!="")
        {
            auto filestream = boost::iostreams::file_source(inputfile, std::ios_base::in | std::ios_base::binary);

            if(!filestream.is_open())
            {
                std::cout << "Could not open " << inputfile << '\n';

                exit(1);
            }

            input.push(filestream);
        }
	}


	void CreateOutputFilter(const std::string& outputfile, boost::iostreams::filtering_ostream& output)
	{
		if(outputfile.substr(outputfile.size()-3) == ".gz")
		{	
			output.push(boost::iostreams::gzip_compressor());
		}

        if(outputfile!="")
        {
            auto filestream = boost::iostreams::file_sink(outputfile, std::ios_base::out | std::ios_base::binary);

            if(!filestream.is_open())
            {
                std::cout << "Could not open " << outputfile << '\n';

                exit(1);
            }

            output.push(filestream);
        }
	}
}


class zifstream: public boost::iostreams::filtering_istream
{
public:
    zifstream(const char* inputfile)
	{
		CreateInputFilter(inputfile, *this);
	}

    zifstream(const std::string& inputfile="")
	{
		CreateInputFilter(inputfile, *this);
	}

    void open(const std::string& inputfile)
    {
        CreateInputFilter(inputfile, *this);
    }
};


class zofstream: public boost::iostreams::filtering_ostream
{
public:
    zofstream(const char* outputfile)
	{
        CreateOutputFilter(outputfile, *this);
	}

    zofstream(const std::string& outputfile)
	{
        CreateOutputFilter(outputfile, *this);
	}

    void open(const std::string& outputfile)
    {
        CreateOutputFilter(outputfile, *this);
    }

    zofstream(){}
};


namespace
{
    int CheckIJformat(const std::string& inputfile)
    {
        ZINPUT(input, inputfile);

        SplitString linedata;

        input >> linedata;

        if(linedata.size() == 6)
        {
            return 0;
        }
        else
        {
            return 1;
        }
    }
}
