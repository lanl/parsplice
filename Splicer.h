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


#ifndef __ParSplice__Splicer__
#define __ParSplice__Splicer__

#include <stdio.h>
#include <limits>
#include <memory>
#include <atomic>
#include <future>
#include <deque>
#include <chrono>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

#include <mpi.h>
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
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/optional.hpp>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <boost/bimap/unordered_multiset_of.hpp>
#include <boost/bimap/vector_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/support/lambda.hpp>
#include <boost/timer/timer.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/random/random_device.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>


#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>



#include "ParSpliceCommon.h"
#include "Validator.h"
#include "ParSpliceWorker.h"
#include "npdds.h"
#include "ldds.h"

#define NCOMMBUFF 16


typedef  boost::bimaps::bimap<
boost::bimaps::unordered_set_of<parsplice::Label>,
boost::bimaps::vector_of< unsigned int >
> TransitionCountsMap;


struct MCNode{
    friend class boost::serialization::access;
    
    TransitionCountsMap counts;
    unsigned int totalCount;
    unsigned int transCount;
    
    MCNode(){
        (*this).totalCount=0;
        (*this).transCount=0;
    };
    
    void incrementTrans(parsplice::Label &l){
        (*this).counts.left[l]++;
        (*this).transCount++;
    };
    
    void incrementTime(unsigned duration){
        (*this).totalCount+=duration;
    };
    
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & totalCount;
        ar & transCount;
        ar & counts;
    }
    
};



struct MC{
   friend class boost::serialization::access;
    
    MC(){
        boost::random::random_device rng;
        rand.seed(rng());
    };
    
    void initFromXML(boost::property_tree::ptree &config){
        forever=config.get<unsigned>("ParSplice.MC.Forever",1000);
    };
    
    void incrementTrans(parsplice::Label &lb, parsplice::Label &le){
        nodes[lb].incrementTrans(le);
    };
    
    void incrementTime(parsplice::Label &lb, unsigned duration){
        nodes[lb].incrementTime(duration);
    };
    
    void update(parsplice::Trajectory &t){
        bool first=true;
        unsigned previous=0;
        
        //std::cout<<"-----"<<std::endl;
        
        
        /*
        //One trajectory contributes to only one state. k_{i->j}=k_{i->*} * branching ratio
        
        //begin
        auto itb=t.visits.begin();
        parsplice::Label lb=itb->label;
        //end
        auto ite=t.visits.rbegin();
        parsplice::Label le=ite->label;
        
        //increment time in the initial state
        unsigned duration=0;
        for(std::list<parsplice::Visit>::iterator it=t.visits.begin();it!=t.visits.end();it++){
            if(it->label==lb){
                duration+=it->duration;
            }
        };
        //unsigned duration=itb->duration;
        incrementTime(lb,duration);
        
        //increment count to the final state
        if(lb!=le){
            incrementTrans(lb,le);
        }
        */
        
        //cannot use this except if we save every intermediate state or if we use coarse trajectories in amdf
        for(std::list<parsplice::Visit>::iterator it=t.visits.begin();it!=t.visits.end();it++){
            //std::cout<<it->label<<" "<<it->duration<<std::endl;
            parsplice::Label l=it->label;
            unsigned duration=it->duration;
            nodes[l].incrementTime(duration);
            if(l==1){
                //std::cout<<"update "<<l<<" "<<duration<<" "<<nodes[l].totalCount<<std::endl;
            }
            //std::cout<<previous<<std::endl;
            if(!first){
                nodes[previous].incrementTrans(l);
                if(previous==1){
                    //std::cout<<"update t "<<previous<<" "<<l<<" "<<nodes[previous].transCount<<" "<<double(nodes[previous].totalCount)/nodes[previous].transCount<<std::endl;
                }
            }
            previous=l;
            first=false;
        };
        
        
    };
    
    //this runs in about (3e-7 + 1e-9*Npath) sec.
    inline std::pair<parsplice::Label,unsigned> sampleBKL(parsplice::Label &lb, double &w){
        
        double r1=uniform(rand);
        double r2=2.;
        
        //this can be set to an arbitarty "big" number
        unsigned duration=forever;
        parsplice::Label final=lb;
        //unknown state, or state without any known exit path
        if( nodes.count(lb)==0 or nodes[lb].transCount==0){
            return std::make_pair(lb,duration);
        }
        else{
        
            MCNode &n=nodes[lb];
            
            //std::cout<<" TC "<<n.transCount<<std::endl;
            double rN=r1*n.transCount;
            double rAccum=0.;
            for(TransitionCountsMap::right_map::const_iterator it=n.counts.right.begin();it!=n.counts.right.end();it++){
                rAccum+=it->first;
                //std::cout<<"    "<<rAccum<<" "<<rN<<std::endl;
                if(rAccum>=rN){
                    final=it->second;
                    break;
                }
            }
            //update the weigth of the trajectory
            w*=double(n.transCount-0.5)/n.transCount;
            
            duration=r2*unsigned(n.totalCount/double(n.transCount));
            
            //std::cout<<lb<<"    "<<n.totalCount<<" "<<n.transCount<<std::endl;
            
            return std::make_pair(final,duration);
        }
    };
    
    //this runs in about (3e-7 + 1e-9*Npath) sec.
    inline std::pair<parsplice::Label,unsigned> sampleBKLV(parsplice::Label &lb, double &w){
        
        double r1=uniform(rand);
        double r2=2.;
        
        //this can be set to an arbitarty "big" number
        unsigned duration=forever;
        parsplice::Label final=lb;
        //unknown state, or state without any known exit path
        if( nodes.count(lb)==0 or nodes[lb].transCount==0){
            return std::make_pair(lb,duration);
        }
        else{
            
            MCNode &n=nodes[lb];
            BOOST_LOG_SEV(lg,info)<<lb<<"    "<<n.totalCount<<" "<<n.transCount<<std::endl;
            
            //std::cout<<" TC "<<n.transCount<<std::endl;
            double rN=r1*n.transCount;
            double rAccum=0.;
            for(TransitionCountsMap::right_map::const_iterator it=n.counts.right.begin();it!=n.counts.right.end();it++){
                rAccum+=it->first;
                BOOST_LOG_SEV(lg,info) <<rAccum<<" "<<rN;
                
                //std::cout<<"    "<<rAccum<<" "<<rN<<std::endl;
                if(rAccum>=rN){
                    final=it->second;
                    break;
                }
            }
            //update the weigth of the trajectory
            w*=double(n.transCount-0.5)/n.transCount;
            
            duration=r2*unsigned(n.totalCount/double(n.transCount));
            
            BOOST_LOG_SEV(lg,info)<<"DURATION    "<<duration<<" "<<w<<std::endl;
            
            //std::cout<<lb<<"    "<<n.totalCount<<" "<<n.transCount<<std::endl;
            
            return std::make_pair(final,duration);
        }
    };
    
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & forever;
        ar & nodes;
    }
    
    
    boost::unordered_map< parsplice::Label, MCNode > nodes;
    boost::random::mt11213b rand;
    boost::random::uniform_01<> uniform;
    boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;
    unsigned forever;
};





class ParSpliceSplicer{
    friend class boost::serialization::access;
    
public:
    
    
    ParSpliceSplicer(MPI_Comm splicerComm_, MPI_Comm dbComm_, std::set<int> children_ , boost::property_tree::ptree &config, bool threadMultiple_){
        BOOST_LOG_SEV(lg,debug) <<"#Splicer() ";
        splicerComm=splicerComm_;
        dbComm=dbComm_;
        threadMultiple=threadMultiple_;
        usedBlocks=unusedBlocks=wastedBlocks=0;
        carryOverTime=0;
        
        children=children_;
        globalValidator.init(children,1);
        
        start=std::chrono::high_resolution_clock::now();
        
        
        bool restartFromCheckPoint=false;
        //check if checkpoint files are present
        restartFromCheckPoint= boost::filesystem::exists("./traj.out.chk") && boost::filesystem::exists("./times.out.chk") && boost::filesystem::exists("./Splicer.chk");
        
        BOOST_LOG_SEV(lg,debug) <<"#Splicer() restartFromCheckPoint: "<<restartFromCheckPoint;
        
        if(restartFromCheckPoint){
            boost::filesystem::copy_file("./traj.out.chk","./traj.out",boost::filesystem::copy_option::overwrite_if_exists);
            boost::filesystem::copy_file("./times.out.chk","./times.out",boost::filesystem::copy_option::overwrite_if_exists);
            outTraj.open("./traj.out", std::ios::app);
            outTime.open("./times.out", std::ios::app);
        }
        else{
            outTraj.open("./traj.out", std::ios::out);
            outTime.open("./times.out", std::ios::out);
        }
        
        
        
        
        
        //initialize from config data
        predictionTime=config.get<unsigned>("ParSplice.MC.PredictionTime",1000);
        wMin=config.get<double>("ParSplice.MC.wMin",0.1);
        
        broadcastDelay=std::chrono::milliseconds( config.get<unsigned>("ParSplice.Splicer.BroadcastDelay",100) );
        reportDelay=std::chrono::milliseconds( config.get<unsigned>("ParSplice.Splicer.ReportDelay",100) );
        checkpointDelay=std::chrono::milliseconds( config.get<unsigned>("ParSplice.Splicer.CheckpointDelay",100) );
        
        
        BOOST_LOG_SEV(lg,debug) <<"#Splicer(): Creating databases ";
        
        //create the databases
        std::map<unsigned int, bool> dbType;
        dbType[1]=true;
        
        std::string sdbType=config.get<std::string>("ParSplice.Splicer.DB.StateDB.DBType");
        std::string base=config.get<std::string>("ParSplice.Splicer.DB.StateDB.BaseName");
        std::string home=config.get<std::string>("ParSplice.Splicer.DB.StateDB.Home");
        boost::trim(sdbType);
        boost::trim(base);
        boost::trim(home);
        unsigned long footprint=config.get<unsigned long>("ParSplice.Splicer.DB.StateDB.Footprint");
        std::string rs=boost::str(boost::format("%1%" ) % 0 );
        boost::replace_all(home, "__RANK__", rs);
        int parent=-1;
        
        
        if(sdbType.find("BDB")!=std::string::npos){
            BOOST_LOG_SEV(lg,trace) <<"#Splicer sdb is DBD";
            sdb=new dds<BDBLocalDataStore>(dbComm,parent,home,base,footprint,dbType,threadMultiple);
        }
        if(sdbType.find("STL")!=std::string::npos){
            BOOST_LOG_SEV(lg,trace) <<"#Splicer sdb is STL";
            sdb=new dds<STLLocalDataStore>(dbComm,parent,home,base,footprint,dbType,threadMultiple);
        }
        
        markovChainEstimator.initFromXML(config);
        
       
        if(restartFromCheckPoint){
            BOOST_LOG_SEV(lg,debug) <<"#Splicer(): Restoring checkpoing ";
            
            // create and open an archive for input
            std::ifstream ifs("Splicer.chk");
            boost::archive::text_iarchive ia(ifs);
            // read class state from archive
            ia >> *this;
            // archive and stream closed when destructors are called
        }
        else{
            //initialize the db and trajectory
            
            parsplice::Label label;
            parsplice::Rd ref;
            driver.ref(label, ref);
            
            sdb->put(DDS_PUT | DDS_FORWARD, 1, label, ref);
            
            parsplice::Visit v;
            v.label=label;
            v.duration=0;
            officialTrajectory.appendVisit(v);
            
            
        }
    };
    
    
    
    
    void server(){
        BOOST_LOG_SEV(lg,debug) <<"#Splicer::server ";
        
        
        //std::vector<parsplice::Label> uintCommBuff(BUFFER_SIZE);
        std::vector< std::vector<parsplice::Label> > uintCommBuffv(NCOMMBUFF);
        for(int i=0;i<NCOMMBUFF;i++){
            uintCommBuffv[i]=std::vector<parsplice::Label>(BUFFER_SIZE);
        }
        
        
        std::vector<parsplice::Label> uintCommBuff2(BUFFER_SIZE);
        MPI_Status status;
        MPI_Request req;
        int completed=0;
        
        
        std::chrono::high_resolution_clock::time_point lastBroadcast=std::chrono::high_resolution_clock::now();
        std::chrono::high_resolution_clock::time_point lastReport=std::chrono::high_resolution_clock::now();
        std::chrono::high_resolution_clock::time_point lastCheckpoint=std::chrono::high_resolution_clock::now();
        std::map<parsplice::Label, unsigned > predictions;
        

        std::map<parsplice::Label, std::multiset<unsigned> > predm;
        
        //post a number of non-blocking receives
        MPI_Request reqv[NCOMMBUFF];
        for(int i=0;i<NCOMMBUFF;i++){
            MPI_Irecv(&(uintCommBuffv[i][0]),BUFFER_SIZE,MPI::UNSIGNED,MPI::ANY_SOURCE,MPI_ANY_TAG,splicerComm,&reqv[i]);
        }
        
        //post an initial non-blocking receive
        //MPI_Irecv(&(uintCommBuff[0]),BUFFER_SIZE,MPI::UNSIGNED,MPI::ANY_SOURCE,MPI_ANY_TAG,splicerComm,&req);
        
        while(true){
            int nMessagesProcessed=0;
            
            for(int i=0;i<NBUFFERSERV;i++){
                completed=0;
                
                while( MPI_Test(&reqv[i],&completed,&status), completed){
                    //test for messages from the outside
                    if(completed){
                        nMessagesProcessed++;
                        //process message
                        int peer=status.MPI_SOURCE;
                        int tag=status.MPI_TAG;
                        int count=0;
                        BOOST_LOG_SEV(lg,trace) <<"#Splicer::server Received message "<<tag<<" from "<<peer<<" in buffer "<<i;
                        MPI_Get_count( &status, MPI::UNSIGNED, &count );
                        
                        
                        
                        if(tag == PARSPLICE_SEGMENTS){
                            BOOST_LOG_SEV(lg,trace) <<"#Splicer::server PARSPLICE_SEGMENTS ";
                            //this is a completed batch of segments
                            int ii=0;
                            parsplice::TrajDB batch;
                            batch=deserializeDB(uintCommBuffv[i], ii);
                            logTrajDBTr(lg,batch);
                            globalValidator.validate(peer,batch);
                            parsplice::TrajDB ret=globalValidator.release();
                            merge(validatedDB,ret);
                        }
                        //post the next receive
                        completed=0;
                        MPI_Irecv(&(uintCommBuffv[i][0]),BUFFER_SIZE,MPI::UNSIGNED,MPI::ANY_SOURCE,MPI_ANY_TAG,splicerComm,&reqv[i]);
                        //this will process receives with higher priority
                    }
                }
            }
            while(sdb->serverSingle()){};
            
            //update the database and Markov chain
            update();
            //splice official trajectory
            bool newEnd=splice();
            //clear up the prediction if the trajectory move to a new state
            if(newEnd){
                //predictions.clear();
                predm.clear();
            }
            //generate a virtual trajectory
            predictions.clear();
            KMCSchedule(predictions);
            for(auto it=predictions.begin();it!=predictions.end();it++){
                predm[it->first].insert(it->second);
            }
            
            
           
            //buffer the predictions and broadcast only if delay has passed since the last time
            if(std::chrono::high_resolution_clock::now() - lastBroadcast> broadcastDelay  ){
                BOOST_LOG_SEV(lg,debug) <<"#Splicer::server broadcasting ";
        
                /*
                std::vector<unsigned> pred;
                
                for(auto it=predictions.begin();it!=predictions.end();it++){
                    pred.push_back(it->first);
                    pred.push_back(it->second);
                    BOOST_LOG_SEV(lg,trace) <<it->first<<" "<<it->second;
                }
                //send to children
                for(auto it=children.begin();it!=children.end();it++){
                    MPI_Send(&(pred[0]),pred.size(),MPI::UNSIGNED,*it,PARSPLICE_TASKS,splicerComm);
                }
                
                predictions.clear();
                */
                
                
                //greedy approach...make this faster
                std::map<parsplice::Label, unsigned > counts;
                int nKeep=0;
                
                for(auto it=predm.begin();it!=predm.end();it++){
                    //counts[it->first]=0;
                    it->second.erase(0);
                }
                
                int granularity=60;
                while(nKeep<predictionTime){
                    unsigned max=0;
                    for(auto it=predm.begin();it!=predm.end();it++){
                        max=std::max(max,unsigned(it->second.size()));
                    }
                    for(auto it=predm.begin();it!=predm.end();it++){
                        if(it->second.size()==max){
                            for(int i=0;i<granularity;i++){
                                nKeep++;
                                counts[it->first]+=1;
                                it->second.erase(counts[it->first]);
                            }
                        }
                    }
                }
                
                std::vector<unsigned> pred;
                for(auto it=counts.begin();it!=counts.end();it++){
                    pred.push_back(it->first);
                    pred.push_back(it->second);
                    BOOST_LOG_SEV(lg,trace) <<it->first<<" "<<it->second;
                }
                //send to children
                for(auto it=children.begin();it!=children.end();it++){
                    
                    //start a non-blocking send
                    MPI_Request r;
                    MPI_Status ss;
                    MPI_Issend(&(pred[0]),pred.size(),MPI::UNSIGNED,*it,PARSPLICE_TASKS,splicerComm, &r);
                    
                    std::chrono::high_resolution_clock::time_point start=std::chrono::high_resolution_clock::now();
                    std::chrono::milliseconds timeout=std::chrono::milliseconds(100);
                    int completed=0;
                    int cancelled=0;
                    
                    while(true){
                        if( std::chrono::high_resolution_clock::now() - start > timeout ){
                            MPI_Cancel(&r);
                            MPI_Wait(&r, &ss);
                            MPI_Test_cancelled(&ss, &cancelled);
                            if(!cancelled){
                                completed=1;
                            }
                            if(cancelled){
                                completed=0;
                            }
                            BOOST_LOG_SEV(lg,error) <<"#Splicer::server cancelling tasks forwarding "<<*it;
                        }
                        else{
                            MPI_Test(&r, &completed, &ss);
                        }
                        
                        
                        if(completed){
                            BOOST_LOG_SEV(lg,debug) <<"#Splicer::server tasks forwarded "<<*it;
                            break;
                        }
                        if(cancelled){
                            BOOST_LOG_SEV(lg,error) <<"#Splicer::server tasks forwarding failed "<<*it;
                            break;
                        }
                        
                    }
                    
                    
                    BOOST_LOG_SEV(lg,debug) <<"#Splicer::server broadcasting done ";
                    //MPI_Send(&(pred[0]),pred.size(),MPI::UNSIGNED,*it,PARSPLICE_TASKS,splicerComm);
                }
                
                predm.clear();
                
                
                
                
                
                lastBroadcast=std::chrono::high_resolution_clock::now();
                
                /*
                KMCScheduleV(predictions);
                BOOST_LOG_SEV(lg,info) <<"SINGLE";
                for(auto it=predictions.begin();it!=predictions.end();it++){
                    BOOST_LOG_SEV(lg,info) <<it->first<<" "<<it->second;
                }
                predictions.clear();
                */
                
            }
            else{
                //BOOST_LOG_SEV(lg,info) <<"#Splicer::server buffering ";
            }
            
            
            //report at given intervals
            if(std::chrono::high_resolution_clock::now() - lastReport> reportDelay  ){
                BOOST_LOG_SEV(lg,debug) <<"#Splicer::server reporting ";
                report();
                output();
                lastReport=std::chrono::high_resolution_clock::now();
                BOOST_LOG_SEV(lg,debug) <<"#Splicer::server reporting done ";
            }
            
            
            //checkpoint at given intervals
            if(std::chrono::high_resolution_clock::now() - lastCheckpoint> checkpointDelay  ){
                BOOST_LOG_SEV(lg,debug) <<"#Splicer::server checkpointing ";
                checkpoint();
                lastCheckpoint=std::chrono::high_resolution_clock::now();
                BOOST_LOG_SEV(lg,debug) <<"#Splicer::server checkpointing done";
            }
        }

        
        
    };

    
    
   
    
    
    
    void KMCSchedule( std::map<parsplice::Label, unsigned > & visits ){
        BOOST_LOG_SEV(lg,trace) <<"#Splicer KMC ";
        
        
        //keep track of which segment in the database we already used
        boost::unordered_map< parsplice::Label, int> consumedSegments;
        
        //start at the current end
        parsplice::Label c=officialTrajectory.back().label;
        double w=1;
        int length=0;
        std::pair<parsplice::Label, unsigned> s;
        {
            //boost::timer::auto_cpu_timer t(3, "KMC scheduling loop: %w seconds\n");
            
            while(true){
                //"virtually" consume segment from the database
                if(segmentDB.count(c)>0 and  consumedSegments[c]< segmentDB[c].size() ){
                    parsplice::Trajectory &t=segmentDB[c][consumedSegments[c]];
                    consumedSegments[c]+=1;
                    BOOST_LOG_SEV(lg,trace) <<"#Splicer KMC virtual consume "<<c<<" "<<t.back().label;
                    c=t.back().label;
                }
                //generate segment with BKL
                else{
                    s=markovChainEstimator.sampleBKL(c,w);
                    BOOST_LOG_SEV(lg,trace) <<"#Splicer KMC generate "<<c<<" "<<s.first<<" "<<s.second<<" "<<w;
                    int duration=s.second;
                    parsplice::Visit v;
                    v.label=c;
                    //do not exceed the prediction horizon
                    v.duration=( duration+length>predictionTime ? predictionTime-length : duration ) ;
                    length+=v.duration;
                    //visits.push_back(v);
                    visits[v.label]+=v.duration;
                    c=s.first;
                    //break if weight becomes too low or we reached the prediction time
                    BOOST_LOG_SEV(lg,trace) <<"#Splicer KMC length "<<length<<" "<<v.duration<<" "<<s.second;
                    if(length>=predictionTime or w<wMin){
                        break;
                    }
                }
            }
        }
        
        
        
        /*
        {
            logging::record rec = lg.open_record();
            if (rec)
            {
                logging::record_ostream strm(rec);
                strm<<"========== PREDICTIONS ==========\n";
                for(auto it=visits.begin();it!=visits.end();it++){
                    strm<<"State: "<<it->label<<" duration "<<it->duration<<"\n";
                }
                strm.flush();
                lg.push_record(boost::move(rec));
            }
        }
        */
    };
    
    
    
    void KMCScheduleV( std::map<parsplice::Label, unsigned > & visits ){
        BOOST_LOG_SEV(lg,info) <<"#Splicer KMCV ";
        
        
        //keep track of which segment in the database we already used
        boost::unordered_map< parsplice::Label, int> consumedSegments;
        
        //start at the current end
        parsplice::Label c=officialTrajectory.back().label;
        double w=1;
        int length=0;
        std::pair<parsplice::Label, unsigned> s;
        {
            //boost::timer::auto_cpu_timer t(3, "KMC scheduling loop: %w seconds\n");
            
            while(true){
                //"virtually" consume segment from the database
                if(segmentDB.count(c)>0 and  consumedSegments[c]< segmentDB[c].size() ){
                    
                    parsplice::Trajectory &t=segmentDB[c][consumedSegments[c]];
                    consumedSegments[c]+=1;
                    BOOST_LOG_SEV(lg,info) <<"#Splicer KMC virtual consume "<<c<<" "<<t.back().label;
                    c=t.back().label;
                                    }
                //generate segment with BKL
                else{
                    s=markovChainEstimator.sampleBKLV(c,w);
                    BOOST_LOG_SEV(lg,info) <<"#Splicer KMC generate "<<c<<" "<<s.first<<" "<<s.second<<" "<<w;
                    int duration=s.second;
                    parsplice::Visit v;
                    v.label=c;
                    //do not exceed the prediction horizon
                    v.duration=( duration+length>predictionTime ? predictionTime-length : duration ) ;
                    length+=v.duration;
                    //visits.push_back(v);
                    visits[v.label]+=v.duration;
                    c=s.first;
                    //break if weight becomes too low or we reached the prediction time
                    BOOST_LOG_SEV(lg,info) <<"#Splicer KMC length "<<length<<" "<<v.duration<<" "<<s.second;
                    if(length>=predictionTime or w<wMin){
                        break;
                    }
                }
            }
        }
        
        
        
        /*
         {
         logging::record rec = lg.open_record();
         if (rec)
         {
         logging::record_ostream strm(rec);
         strm<<"========== PREDICTIONS ==========\n";
         for(auto it=visits.begin();it!=visits.end();it++){
         strm<<"State: "<<it->label<<" duration "<<it->duration<<"\n";
         }
         strm.flush();
         lg.push_record(boost::move(rec));
         }
         }
         */
    };

    
    
    void update(){
        BOOST_LOG_SEV(lg,trace) <<"#Splicer update ";
        //logTrajDB(lg,validatedDB);
        
        if(validatedDB.size()>0){
            
            //update the markovChainEstimator
            for(parsplice::TrajDB::iterator it=validatedDB.begin();it!=validatedDB.end();it++){
                //loop over trajectories
                for(std::deque<parsplice::Trajectory>::iterator itt=it->second.begin();itt!=it->second.end();itt++){
                    parsplice::Trajectory &t=*itt;
                    unusedBlocks+=t.length;
                    markovChainEstimator.update(t);
                }
            }
            //merge new segments into the main database
            merge(segmentDB, validatedDB);
            
            validatedDB.clear();
        }
        
        //logTrajDB(lg,segmentDB);
    };
     
    bool splice(){
        BOOST_LOG_SEV(lg,trace) <<"#Splicer::splice ";
        //add segments to the official trajectory
        parsplice::Label end=officialTrajectory.back().label;
        //BOOST_LOG_SEV(lg,info) <<"#Splicer Splice Current end: "<<end;
        bool endChanged=false;
        parsplice::Label originalEnd=end;
        
        while( segmentDB[end].size()> 0 ){
            
            parsplice::Trajectory t=segmentDB[end].front();
            segmentDB[end].pop_front();
            usedBlocks+=t.length;
            unusedBlocks-=t.length;
            officialTrajectory.splice(t);
            end=officialTrajectory.back().label;
            BOOST_LOG_SEV(lg,trace) <<"#Splicer Splice Current end: "<<end<<" "<<t.length;
        }
        
        return originalEnd !=officialTrajectory.back().label ;
        
    };
    
    void report(){
        BOOST_LOG_SEV(lg,info) <<"#====================================== ";
        BOOST_LOG_SEV(lg,info) <<"#Splicer::report "<<std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now()-start).count();
        parsplice::Label end=officialTrajectory.back().label;
        BOOST_LOG_SEV(lg,info) <<"#Splicer Splice Current end: "<<end;
        
        /*
        BOOST_LOG_SEV(lg,info) <<"#Splicer Database: ";
        logTrajDB(lg,segmentDB);
        BOOST_LOG_SEV(lg,info) <<"#====================================== ";
        */
    };
    
    void output(){
        
        //output timings
        outTime<<std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now()-start).count()+carryOverTime<<" "<<usedBlocks<<" "<<unusedBlocks<<std::endl;
        
        
        //output trajectory
        unsigned int nBlocks=0;
        int iVisit=0;
        parsplice::Visit back=officialTrajectory.back();
        officialTrajectory.visits.pop_back();
        
        for(auto it=officialTrajectory.visits.begin();it!=officialTrajectory.visits.end();it++){
            nBlocks+=it->duration;
            outTraj<<" "<<it->label<<" "<<it->duration<<"\n";
            iVisit++;
        }
        outTraj.flush();
        
        officialTrajectory.clear();
        officialTrajectory.appendVisit(back);
        
    };
    
    void checkpoint(){
        sdb->sync();
        
        outTraj.close();
        outTime.close();
        
        //copy consistent progress files
        boost::filesystem::copy_file("./traj.out","./traj.out.chk",boost::filesystem::copy_option::overwrite_if_exists);
        boost::filesystem::copy_file("./times.out","./times.out.chk",boost::filesystem::copy_option::overwrite_if_exists);
        
        outTraj.open("./traj.out", std::ios::app | std::ios::out);
        outTime.open("./times.out", std::ios::app | std::ios::out);
        
        
        std::ofstream ofs("./Splicer.chk");
        // save data to archive
        {
            boost::archive::text_oarchive oa(ofs);
            // write class instance to archive
            oa << *this;
            // archive and stream closed when destructors are called
        }

        
    };
    
    
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const{
        // note, version is always the latest when saving
        ar  & usedBlocks;
        ar & unusedBlocks;
        long int co=std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now()-start).count()+carryOverTime;
        ar  & co;
        ar & segmentDB;
        ar & markovChainEstimator;
        ar & officialTrajectory;
        
    };
    
    template<class Archive>
    void load(Archive & ar, const unsigned int version){
        ar  & usedBlocks;
        ar & unusedBlocks;
        ar & carryOverTime;
        ar & segmentDB;
        ar & markovChainEstimator;
        ar & officialTrajectory;
    };
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    
    parsplice::TrajDB validatedDB;
    parsplice::Trajectory officialTrajectory;
    
protected:
    parsplice::TrajDB segmentDB;
    std::set<int> children;
    Validator globalValidator;
    
    
    MC markovChainEstimator;
    unsigned long usedBlocks;
    unsigned long unusedBlocks;
    unsigned long wastedBlocks;

    
    int predictionTime;
    double wMin;
    MPI_Comm splicerComm;
    MPI_Comm dbComm;
    
    
    boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;
    
    int carryOverTime;
    
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::milliseconds broadcastDelay;
    std::chrono::milliseconds reportDelay;
    std::chrono::milliseconds checkpointDelay;

    std::ofstream outTraj;
    std::ofstream outTime;
    
    
    AbstractDDS *sdb;

    bool threadMultiple;
    ParSpliceWorkerDriver driver;
};


#endif /* defined(__ParSplice__Splicer__) */
