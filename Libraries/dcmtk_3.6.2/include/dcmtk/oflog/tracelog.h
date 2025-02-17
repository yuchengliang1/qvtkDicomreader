// -*- C++ -*-
// Module:  Log4CPLUS
// File:    tracelog.h
// Created: 1/2009
// Author:  Vaclav Haisman
//
//
// Copyright 2009-2010 Tad E. Smith
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

/** @file */

#ifndef DCMTK_LOG4CPLUS_TRACELOGGER_H
#define DCMTK_LOG4CPLUS_TRACELOGGER_H

#include "dcmtk/oflog/config.h"

#if defined (DCMTK_LOG4CPLUS_HAVE_PRAGMA_ONCE)
#pragma once
#endif

#include "dcmtk/oflog/logger.h"


namespace dcmtk
{
	namespace log4cplus
	{


		/**
		* This class is used to produce "Trace" logging.  When an instance of
		* this class is created, it will log a <code>"ENTER: " + msg</code>
		* log message if TRACE_LOG_LEVEL is enabled for <code>logger</code>.
		* When an instance of this class is destroyed, it will log a
		* <code>"ENTER: " + msg</code> log message if TRACE_LOG_LEVEL is enabled
		* for <code>logger</code>.
		* <p>
		* @see DCMTK_LOG4CPLUS_TRACE
		*/
		class TraceLogger
		{
		public:
			TraceLogger(const Logger& l, const log4cplus::tstring& _msg,
				const char* _file = NULL, int _line = -1)
				: logger(l), msg(_msg), file(_file), line(_line)
			{
				if (logger.isEnabledFor(TRACE_LOG_LEVEL))
					;
					//logger.forcedLog(TRACE_LOG_LEVEL, DCMTK_LOG4CPLUS_TEXT("ENTER: ") + msg, file, line);
					/*
					 * 此处的bug是由于字符集问题导致
					 */
			}

			~TraceLogger()
			{
				if (logger.isEnabledFor(TRACE_LOG_LEVEL))
					;
					//logger.forcedLog(TRACE_LOG_LEVEL, DCMTK_LOG4CPLUS_TEXT("EXIT:  ") + msg, file, line);
					/*
					* 此处的bug是由于字符集问题导致
					*/
			}

		private:
			TraceLogger(TraceLogger const &);
			TraceLogger & operator = (TraceLogger const &);

			Logger logger;
			log4cplus::tstring msg;
			const char* file;
			int line;
		};


	} // log4cplus
} // end namespace dcmtk


#endif // DCMTK_LOG4CPLUS_TRACELOGGER_H
