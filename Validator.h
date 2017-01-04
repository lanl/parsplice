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


#ifndef __ParSplice__Validator__
#define __ParSplice__Validator__

#include <stdio.h>

#include <set>
#include <vector>
#include <list> 
#include <limits>
#include <pthread.h>

#include <mpi.h>
#include "ParSpliceCommon.h"





class Validator{
public:
    
    /*
    Validator(std::set<int> clients_,int batchSize_){
        clients=clients_;
        batchSize=batchSize_;
        ib=0;
        wastedBlocks=0;
        
        
        for(auto it=clients.begin();it!=clients.end();it++){
            Placeholder p;
            p.producerID=*it;
            unvalidated.push_back(p);
        }
        
    };
    */
    
    void init(std::set<int> clients_,int batchSize_){
        clients=clients_;
        batchSize=batchSize_;
        
        ib=0;
        wastedBlocks=0;
        for(auto it=clients.begin();it!=clients.end();it++){
            Placeholder p;
            p.producerID=*it;
            unvalidated.push_back(p);
        }

    }
    void validate(int peer, parsplice::TrajDB &batch){
        BOOST_LOG_SEV(lg,debug) <<"#Validator::validate "<<peer<<" "<<unvalidated.size();
        
        bool phFound=false;
        //for(auto it=unvalidated.rbegin();it!=unvalidated.rend();it++)
        for(auto it=unvalidated.begin();it!=unvalidated.end();it++)
        {
            if(it->producerID==peer and it->db.size()==0){
                it->db=batch;
                phFound=true;
                break;
            }
        }
        
        //placeholder not found
        if(not phFound){
            BOOST_LOG_SEV(lg,warning) <<"#Validator::validate placeholder not found";
            //queue at the end
            Placeholder p;
            p.producerID=peer;
            p.db=batch;
            unvalidated.push_back(p);
        }
        
        //insert new placeholder at the end
        {
            Placeholder p;
            p.producerID=peer;
            unvalidated.push_back(p);
        }
    };
    
    parsplice::TrajDB release(){
        //BOOST_LOG_SEV(lg,info) <<"#Validator::release ";
        /*
        int a;
        do{
            a=0;
            //remove placeholders that have expired (and put them back at the end?)
            if(unvalidated.size()>0){
                boost::timer::cpu_times age=unvalidated.front().timer.elapsed();
                auto nanoseconds = boost::chrono::nanoseconds(age.user);
                auto as=boost::chrono::duration_cast<boost::chrono::seconds>(nanoseconds);
                
                //only empty placeholders can expire
                if(unvalidated.front().db.size()==0){
                    a=as.count();
                }
            }
            if(a>expirationTime){
                unvalidated.pop_front();
            }
        }while(a>expirationTime);
        */
        
        //release validated segments
        parsplice::TrajDB ret;
        
        //std::cout<<"b "<<unvalidated.size()<<std::endl;
        
        /*
        int ic=0;
        for(auto it=unvalidated.begin();it!=unvalidated.end();it++){
            std::cout<<it->producerID<<" "<<it->db.size()<<std::endl;
            ic++;
            if(ic>10){
                break;
            }
        }
        */
        
        /*
        {
            logging::record rec = lg.open_record(keywords::severity = info);
            if (rec)
            {
                logging::record_ostream strm(rec);
                strm<<"====================\n";
                int ic=0;
                for(auto it=unvalidated.begin();it!=unvalidated.end();it++){
                    strm<<it->producerID<<" "<<it->db.size()<<"\n";
                    ic++;
                    if(ic>10){
                        break;
                    }
                }
                strm.flush();
                lg.push_record(boost::move(rec));
            }
            
        }
        */
        
        //std::cout<<"fs "<<unvalidated.front().db.size()<<std::endl;
        while(unvalidated.size()>0 and unvalidated.front().db.size()>0){
            merge(ready, unvalidated.front().db);
            unvalidated.pop_front();
            ib++;
            if(ib==batchSize){
                merge(ret, ready);
                ready.clear();
                ib=0;
                break;
            }
        }
        
        //std::cout<<"e "<<unvalidated.size()<<std::endl;
        return ret;
    };
    
    
    
    
private:
    int ib;
    std::set<int> clients;
    int batchSize;
    parsplice::TrajDB validated;
    std::list< Placeholder > unvalidated;
    parsplice::TrajDB ready;
    unsigned wastedBlocks;
    
    boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;
};



//TimedValidator might be broken.... think the logic through again.



class TimedValidator{
public:
    
    void init(std::set<int> clients_,int delay_){
        clients=clients_;
        
        last=std::chrono::high_resolution_clock::now();
        delay=std::chrono::milliseconds( delay_ );
        
        ib=0;
        wastedBlocks=0;
        for(auto it=clients.begin();it!=clients.end();it++){
            Placeholder p;
            p.producerID=*it;
            unvalidated.push_back(p);
        }
        
    }
    void validate(int peer, parsplice::TrajDB &batch){
        BOOST_LOG_SEV(lg,debug) <<"#Validator::validate "<<peer;
        
        bool phFound=false;
        for(auto it=unvalidated.rbegin();it!=unvalidated.rend();it++){
            if(it->producerID==peer and it->db.size()==0){
                it->db=batch;
                phFound=true;
                break;
            }
        }
        
        //placeholder not found
        if(not phFound){
            BOOST_LOG_SEV(lg,warning) <<"#Validator::validate placeholder not found";
            //queue at the end
            Placeholder p;
            p.producerID=peer;
            p.db=batch;
            unvalidated.push_back(p);
        }
        
        //insert new placeholder at the end
        {
            Placeholder p;
            p.producerID=peer;
            unvalidated.push_back(p);
        }
    };
    
    parsplice::TrajDB release(){
        
        parsplice::TrajDB ret;
        
        //std::cout<<"fs "<<unvalidated.front().db.size()<<std::endl;
        while(unvalidated.size()>0 and unvalidated.front().db.size()>0){
            merge(ready, unvalidated.front().db);
            unvalidated.pop_front();
            ib++;
        }
        
        if(std::chrono::high_resolution_clock::now() - last> delay  ){
            merge(ret, ready);
            ready.clear();
            ib=0;
            last=std::chrono::high_resolution_clock::now();
        }
        
        //std::cout<<"e "<<unvalidated.size()<<std::endl;
        return ret;
    };
    
    
    
    
private:
    int ib;
    std::set<int> clients;
    std::chrono::milliseconds delay;
    std::chrono::high_resolution_clock::time_point last;
    parsplice::TrajDB validated;
    std::list< Placeholder > unvalidated;
    parsplice::TrajDB ready;
    unsigned wastedBlocks;
    
    boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;
};




























/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////






















#endif /* defined(__ParSplice__Validator__) */
