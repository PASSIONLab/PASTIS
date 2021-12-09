/**
 * @file
 *	logger.hpp
 *
 * @author
 *	Oguz Selvitopi
 *
 * @date
 *
 * @brief
 *	Simple logging utility
 *
 * @todo
 *
 * @note
 *	
 */

#pragma once

#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#define PASTIS_VERBOSITY 3



namespace
pastis
{


// Can use a file or standard out/err for logging
class
Logger
{
	
public:

	enum class LogLevel {INFO, INFO2, INFO3, WARNING, ERROR};

	// const std::string infostr{"[INFO1]"}, infostrwrap{"\u2593"};
	// const std::string info2str{"[INFO2]"}, info2strwrap{"\u2592"}; 
	// const std::string info3str{"[INFO3]"}, info3strwrap{"\u2591"};
	// const std::string warnstr{"[WARN]"}, warnstrwrap{"\u2212"};
	// const std::string errstr{"[ERROR]"}, errstrwrap{"\u2217"};
	const std::string infostr{"[INFO1]"}, infostrwrap{"INF1"};
	const std::string info2str{"[INFO2]"}, info2strwrap{"INF2"}; 
	const std::string info3str{"[INFO3]"}, info3strwrap{"INF3"};
	const std::string warnstr{"[WARN]"}, warnstrwrap{"WARN"};
	const std::string errstr{"[ERROR]"}, errstrwrap{"ERR"};	


	static
	std::shared_ptr<Logger>
	instantiate(bool active,
				std::ostream &out = std::cout)
	{
		if (!logger_)
			logger_ = std::shared_ptr<Logger>(new Logger(active, out));
		return logger_;
	}
	


	static
	std::shared_ptr<Logger>
	instantiate(bool active,
				const std::string &fname)
	{
		if (!logger_)
			logger_ = std::shared_ptr<Logger>(new Logger(active, fname));
		return logger_;
	}


	
	void
	log (const std::string &s,
		 LogLevel ll = LogLevel::INFO)
	{
		if (out_ == nullptr)
			return;

		// file and line number
		// std::stringstream ss1;
		// ss1 << "(" __FILE__ << ":" << std::setw(4) << __LINE__ << ")";

		// time
		std::stringstream ss2;
		time(&tlog_);
		struct tm *tinfo = localtime(&tlog_);
		ss2 << std::setfill('0') << std::setw(4) << tinfo->tm_year+1900 << "-"
			<< std::setw(2) << tinfo->tm_mon+1 << "-"
			// << std::setw(2) << tinfo->tm_mday << "\u22C6"
			<< std::setw(2) << tinfo->tm_mday << "|"
			<< std::setw(2) << tinfo->tm_hour << ":"
			<< std::setw(2) << tinfo->tm_min << ":"
			<< std::setw(2) << tinfo->tm_sec;		
		
		switch (ll)
		{
			#if PASTIS_VERBOSITY > 0
		case LogLevel::INFO:
			(*out_) << infostr << " " << infostrwrap << " "
					<< ss2.str() << " " << infostrwrap << " "
					// << ss1.str() << " " << infostrwrap << " "
					<< s << std::endl;
			break;
			#if PASTIS_VERBOSITY > 1
		case LogLevel::INFO2:
			(*out_) << info2str << " " << info2strwrap << " "
					<< ss2.str() << " " << info2strwrap << " "
					// << ss1.str() << " " << info2strwrap << " "
					<< s << std::endl;
			break;
			#if PASTIS_VERBOSITY > 2
		case LogLevel::INFO3:
			(*out_) << info3str << " " << info3strwrap << " "
					<< ss2.str() << " " << info3strwrap << " "
					// << ss1.str() << " " << info3strwrap << " "
					<< s << std::endl;
			break;
			#endif
			#endif
			#endif
		case LogLevel::WARNING:
			(*out_) << warnstr << " " << warnstrwrap << " "
					<< ss2.str() << " " << warnstrwrap << " "
					// << ss1.str() << " " << warnstrwrap << " "
					<< s << std::endl;
			break;
		case LogLevel::ERROR:
			(*out_) << errstr << " " << errstrwrap << " "
					<< ss2.str() << " " << errstrwrap << " "
					// << ss1.str() << " " << errstrwrap << " "
					<< s << std::endl;
			break;
		default:
			break;
		}

		return;
	}



	void
	debug (const std::string &s)
	{
		(*dbgf_) << s << std::endl;
	}
	


	void
	set_dbgfile (const std::string &fname)
	{
		if (dbgf_ != nullptr)
			dbgf_->close();
		dbgf_->open(fname, std::ofstream::out);
	}



	void
	close_dbgfile ( )
	{
		if (dbgf_ != nullptr)
			dbgf_->close();
	}



	
	~Logger ()
	{
		if (outf_ != nullptr)
		{
			outf_->close();
			delete outf_;
		}
		if (dbgf_ != nullptr)
		{
			dbgf_->close();
			delete dbgf_;
		}
	}
	


private:

	static std::shared_ptr<Logger>	 logger_;
	std::ostream					*out_;
	std::ofstream					*outf_;
	std::time_t						 tlog_;
	std::ofstream                   *dbgf_;
	
	Logger (bool active,
			std::ostream &out = std::cout)
	{
		out_ = nullptr;
		outf_ = nullptr;
		if (active)
			out_ = &out;
		dbgf_ = new std::ofstream;
	}
	

	Logger (bool active,
			const std::string &fname)
	{
		out_ = nullptr;
		outf_ = nullptr;
		if (active)
		{
			outf_ = new std::ofstream;
			outf_->open(fname);
			out_ = outf_;
		}
		dbgf_ = new std::ofstream;
	}
};


}
