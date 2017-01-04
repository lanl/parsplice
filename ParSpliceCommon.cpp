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


#include "ParSpliceCommon.h"


void imerge(parsplice::TrajDB &source, parsplice::TrajDB &supp){
    
    for(parsplice::TrajDB::iterator it=source.begin();it!=source.end();it++){
        //loop over trajectories
        for(std::deque<parsplice::Trajectory>::iterator itt=it->second.begin();itt!=it->second.end();itt++){
            std::cout<<"=========source===========\n";
            for(std::list<parsplice::Visit>::iterator ittt=itt->visits.begin();ittt!=itt->visits.end();ittt++){
                std::cout<<" -- ( "<<ittt->label<<" ) "<<ittt->duration;
            }
            std::cout<<"\n";
        }
    }
    
    for(parsplice::TrajDB::iterator it=supp.begin();it!=supp.end();it++){
        //loop over trajectories
        for(std::deque<parsplice::Trajectory>::iterator itt=it->second.begin();itt!=it->second.end();itt++){
            std::cout<<"=========supp===========\n";
            for(std::list<parsplice::Visit>::iterator ittt=itt->visits.begin();ittt!=itt->visits.end();ittt++){
                std::cout<<" -- ( "<<ittt->label<<" ) "<<ittt->duration;
            }
            std::cout<<"\n";
        }
    }
    
    
    //loop over initial states
    for(parsplice::TrajDB::iterator it=supp.begin();it!=supp.end();it++){
        std::cout<<"o "<<source.size()<<" "<<supp.size()<<std::endl;
        //loop over trajectories
        for(std::deque<parsplice::Trajectory>::iterator itt=it->second.begin();itt!=it->second.end();itt++){
            std::cout<<"i "<<it->first<<" "<<itt->empty()<<std::endl;
            parsplice::Trajectory &t=*itt;
            
            if(not t.empty()){
                
                parsplice::Label lb=t.front().label;
                std::cout<<"ii "<<it->first<<" label "<<lb<<" "<<source.count(lb)<<std::endl;
                
                //if trajectory exist in that bin, try to splice at back. Otherwise, add
                if(source.count(lb)>0 and not source[lb].empty()){
                    if(not source[lb].back().splice(t)){
                        source[lb].push_back(t);
                    }
                }
                else{
                    std::cout<<"iii "<<it->first<<" label "<<lb<<" "<<source.count(lb)<<std::endl;
                    source[lb]=std::deque<parsplice::Trajectory>();
                    source[lb].push_back(t);
                }
            }
        }
    }
};


void merge(parsplice::TrajDB &source, parsplice::TrajDB &supp){
    
    //loop over initial states
    for(parsplice::TrajDB::iterator it=supp.begin();it!=supp.end();it++){
        //loop over trajectories
        for(std::deque<parsplice::Trajectory>::iterator itt=it->second.begin();itt!=it->second.end();itt++){
            parsplice::Trajectory &t=*itt;
            if(not t.empty()){
                parsplice::Label lb=t.front().label;
                //if trajectory exist in that bin, try to splice at back. Otherwise, add
                if(source.count(lb)>0 and not source[lb].empty()){
                    if(not source[lb].back().splice(t)){
                        source[lb].push_back(t);
                    }
                }
                else{
                    source[lb]=std::deque<parsplice::Trajectory>();
                    source[lb].push_back(t);
                }
            }
        }
    }
};



void serializeDB(parsplice::TrajDB &source, std::vector<unsigned int> &v, int &i){
    
    int i0=i;
    //placeholder for the number of trajectories in the db
    v[i++]=0;
    parsplice::Trajectory t;
    
    
    //loop over initial states
    for(parsplice::TrajDB::iterator it=source.begin();it!=source.end();it++){
        //loop over trajectories
        for(std::deque<parsplice::Trajectory>::iterator itt=it->second.begin();itt!=it->second.end();itt++){
            itt->serialize(v,i);
            v[i0]++;
        }
    }
};


parsplice::TrajDB deserializeDB( std::vector<unsigned int> &v, int &i){
    parsplice::TrajDB db;
    parsplice::Trajectory t;
    int nt=v[i++];
    
    for(int j=0;j<nt;j++){
        t.clear();
        t.deserialize(v, i);
        
        parsplice::Label l=t.front().label;
        db[l].push_back(t);
    }
    
    
    return db;
};


void logTrajDB(boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg, parsplice::TrajDB &db){
    logging::record rec = lg.open_record(keywords::severity = info);
    if (rec)
    {
        logging::record_ostream strm(rec);
        strm<<"========== DB ==========\n";
        
        //loop over initial states
        for(parsplice::TrajDB::iterator it=db.begin();it!=db.end();it++){
            //loop over trajectories
            for(std::deque<parsplice::Trajectory>::iterator itt=it->second.begin();itt!=it->second.end();itt++){
                strm<<"====================\n";
                strm<<itt->length<<"\n";
                for(std::list<parsplice::Visit>::iterator ittt=itt->visits.begin();ittt!=itt->visits.end();ittt++){
                    strm<<" -- ( "<<ittt->label<<" ) "<<ittt->duration;
                }
                strm<<"\n";
            }
        }
        strm.flush();
        lg.push_record(boost::move(rec));
    }
    
};

void logTrajDBTr(boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg, parsplice::TrajDB &db){
    logging::record rec = lg.open_record(keywords::severity = trace);
    if (rec)
    {
        logging::record_ostream strm(rec);
        strm<<"========== DB ==========\n";
        
        //loop over initial states
        for(parsplice::TrajDB::iterator it=db.begin();it!=db.end();it++){
            //loop over trajectories
            for(std::deque<parsplice::Trajectory>::iterator itt=it->second.begin();itt!=it->second.end();itt++){
                strm<<"====================\n";
                strm<<itt->length<<"\n";
                for(std::list<parsplice::Visit>::iterator ittt=itt->visits.begin();ittt!=itt->visits.end();ittt++){
                    strm<<" -- ( "<<ittt->label<<" ) "<<ittt->duration;
                }
                strm<<"\n";
            }
        }
        strm.flush();
        lg.push_record(boost::move(rec));
    }
    
};

