/*
 Copyright (c) 2016, Los Alamos National Security, LLC
 All rights reserved.
 Copyright 2016. Los Alamos National Security, LLC. This software was produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL.
 
 Additionally, redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 1.      Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 2.      Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 3.      Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL, the U.S. Government, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


#ifndef __ParSplice__Common__
#define __ParSplice__Common__

#include <stdio.h>
#include <vector>
#include <list>
#include <deque>
#include <memory>
#include <atomic>
#include <future>


#include <boost/unordered_map.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/utility/setup/settings.hpp>
#include <boost/log/utility/setup/from_settings.hpp>

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/timer/timer.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_01.hpp>

#include <mpi.h>

#include "ParSpliceTypes.h"


#define PARSPLICE_TASKS 1 << 1
#define PARSPLICE_SEGMENTS 1 << 2

namespace logging = boost::log;
namespace sinks = boost::log::sinks;
namespace src = boost::log::sources;
namespace expr = boost::log::expressions;
namespace attrs = boost::log::attributes;
namespace keywords = boost::log::keywords;
using namespace logging::trivial;

#define BUFFER_SIZE 100000000



void imerge(parsplice::TrajDB &source, parsplice::TrajDB &supp);
void merge(parsplice::TrajDB &source, parsplice::TrajDB &supp);

void serializeDB(parsplice::TrajDB &source,  std::vector<unsigned int> &v, int &i);
parsplice::TrajDB deserializeDB(std::vector<unsigned int> &v, int &i);

struct Placeholder{
    boost::timer::cpu_timer timer;
    int producerID;
    parsplice::TrajDB db;
    Placeholder(){
        timer.start();
    };
    
};



void logTrajDB(boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg, parsplice::TrajDB &db);

void logTrajDBTr(boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg, parsplice::TrajDB &db);

#endif /* defined(__ParSplice__Common__) */
