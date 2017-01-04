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


#include <iostream>


#include <boost/timer/timer.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/info_parser.hpp>



#include "Splicer.h"
#include "NodeManager.h"




int main(int argc,  char * argv[]) {
    std::cout << "ParSPLICE\n";
    
   
    
    int provided;
    MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &provided );
    std::cout<<"MPI THREAD LEVEL "<<provided<<" REQUESTED "<<MPI_THREAD_MULTIPLE<<std::endl;
    
    
    parsplice::NoSpinLock mpiLock;
    parsplice::AbstractSpinLock *commLock=static_cast< parsplice::AbstractSpinLock * >(&mpiLock);
    
    
    
    
    int rank=0;
    int nRanks=1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

   
    // Create empty property tree object
    boost::property_tree::ptree tree;
    // Parse the XML into the property tree.
    boost::property_tree::read_xml("./input/ps-config.xml", tree);
    //boost::property_tree::read_info("./input/ps-config.info", tree);
    
    
    bool threadMultiple=tree.get("ParSplice.ThreadMultiple", true);
    if(threadMultiple and (provided!=MPI_THREAD_MULTIPLE)){
        std::cerr<<" ParSPLICE REQUIRES MPI_THREAD_MULTIPLE"<<std::endl;
        return 0;
    }
    
    int arity=tree.get("ParSplice.Topology.Arity", 2);
    
    
    
    std::string prefix="rank.__RANK__";
    std::string rs=boost::str(boost::format("%1%" ) % rank );
    boost::replace_all(prefix, "__RANK__", rs);
    
    int loggingSeverity=tree.get("ParSplice.Logging.Severity", 0);

    
    {
        std::string file=str(boost::format("%s.log") % prefix);
        logging::settings setts;
        
        setts["Core"]["DisableLogging"] = false;
        setts["Sinks.Destination"] = "TextFile";
        setts["Sinks.File.Destination"] = "TextFile";
        setts["Sinks.File.FileName"] = file;
        setts["Sinks.File.AutoFlush"] = true;
        setts["Sinks.File.Format"] = "[ %TimeStamp% ] %Severity%: %Message%";
        
        //keywords::format = "[%TimeStamp%]: %Message%"
        logging::register_simple_formatter_factory< boost::log::trivial::severity_level, char >("Severity");

        logging::add_common_attributes();
        
        
        logging::init_from_settings(setts);

        switch(loggingSeverity){
            case 6:
                logging::core::get()->set_filter
                (
                 logging::trivial::severity >= logging::trivial::fatal
                 );
                break;
            case 5:
                logging::core::get()->set_filter
                (
                 logging::trivial::severity >= logging::trivial::error
                 );
                break;
            case 4:
                logging::core::get()->set_filter
                (
                 logging::trivial::severity >= logging::trivial::warning
                 );
                break;
            case 3:
                logging::core::get()->set_filter
                (
                 logging::trivial::severity >= logging::trivial::info
                 );
                break;
            case 2:
                logging::core::get()->set_filter
                (
                 logging::trivial::severity >= logging::trivial::debug
                 );
                break;
            case 1:
                logging::core::get()->set_filter
                (
                 logging::trivial::severity >= logging::trivial::trace
                 );
                break;
            default:
                //setts["Core"]["DisableLogging"] = true;
                break;
        }
        
    }
    boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;
    
    
    //identify parent and children using a k-ary tree
    std::set<int> children;
    int parent=floor((rank-1.0)/arity);
    
    for(int c=0;c<arity;c++){
        int ic=arity*rank+1+c;
        if(ic<nRanks and ic>=0){
            children.insert(ic);
        }
    }
    std::vector<int> ranks(nRanks,0);
    for(int i=0;i<nRanks;i++){
        ranks[i]=i;
    }
    
    bool iAmRoot=(parent==-1);
    bool iAmLeaf=(children.size()==0);
    int depth=0;
    int s=0;
    while(s<rank){
        depth++;
        s+=pow(arity,depth);
    }
    
    
    /*
    std::cout<<"RANK "<<rank<<" PARENT "<<parent<<" DEPTH "<<depth<<" CHILDREN ";
    for(std::set<int>::iterator it=children.begin();it!=children.end();it++){
        std::cout<<*it<<" ";
    }
    std::cout<<std::endl;
    */
    
    
    
    
    MPI_Comm splicerComm;
    MPI_Comm ddsComm;
    

    //create a validator comm that is just a copy of world
    MPI_Comm_dup(MPI_COMM_WORLD, &splicerComm);
    //create a second dispatch comm that is just a copy of world
    //the dispatchers and producers can talk over that comm
    MPI_Comm_dup(MPI_COMM_WORLD, &ddsComm);

    
    
    
#if 1
    if(depth==0){
        //splicer
        ParSpliceSplicer splicer(splicerComm, ddsComm, children, tree, threadMultiple);
        splicer.server();
    }
    
    
    if(depth>0){
        NodeManager manager(splicerComm, ddsComm, rank,children, parent, tree, threadMultiple);
        manager.server();
        
        while(true){};
    }
#endif
    
    
  
   
    MPI_Finalize();
    return 0;
}





































#if 0

ParSpliceWorkerDriver driver;
parsplice::Label label;
parsplice::Rd ref;
driver.ref(label, ref);

std::cout<<label<<std::endl;

if(depth==0){
    //splicer validator db
    
    
    std::string homeDir=str(boost::format("./r%i/") % rank);
    std::string baseName="test";
    unsigned long scratchSize=10000000;
    //dds<STLLocalDataStore> ds(ddsComm,parent,homeDir,baseName,scratchSize,dbType);
    dds<BDBLocalDataStore> ds(ddsComm,parent,homeDir,baseName,1,dbType);
    
    
    
    
    for(int i=0;i<1000;i++){
        boost::timer::auto_cpu_timer t;
        
        label+=1;
        ds.put(DDS_PUT | DDS_FORWARD, 1, label, ref);
        
        /*
         std::cout<<label<<std::endl;
         parsplice::Rd reft;
         std::shared_future<Rd> rf= ds.get( DDS_GET | DDS_FORWARD, 1, label);
         while( rf.wait_for(std::chrono::seconds(0)) != std::future_status::ready){
         
         }
         
         
         Rd data=rf.get();
         std::cout<<"DONE "<<data.size()<<std::endl;
         */
    }
    
    
    
    
    /*
     std::cout<<label<<std::endl;
     parsplice::Rd reft;
     std::shared_future<Rd> rf= ds.get( DDS_GET | DDS_FORWARD, 1, label);
     
     
     while( rf.wait_for(std::chrono::seconds(0)) != std::future_status::ready){
     
     }
     
     Rd data=rf.get();
     for(int i=0;i<data.size();i++){
     std::cout<<int(data[i]);
     }
     */
    
    
    while(true){
        
    }
    
    
    
    
}


if(depth>0){
    std::string homeDir=str(boost::format("./r%i/") % rank);
    std::string baseName="test";
    unsigned long scratchSize=10000000;
    //dds<STLLocalDataStore> ds(ddsComm,parent,homeDir,baseName,scratchSize,dbType);
    dds<BDBLocalDataStore> ds(ddsComm,parent,homeDir,baseName,1,dbType);
    
    
    
    //std::this_thread::sleep_for(std::chrono::seconds(3));
    
    for(int i=0;i<1000;i++){
        boost::timer::auto_cpu_timer t;
        
        label+=1;
        //ds.put(DDS_PUT | DDS_FORWARD, 1, label, ref);
        
        
        std::cout<<label<<std::endl;
        parsplice::Rd reft;
        std::shared_future<Rd> rf= ds.get( DDS_GET | DDS_FORWARD, 1, label);
        while( rf.wait_for(std::chrono::seconds(0)) != std::future_status::ready){
            
        }
        
        
        Rd data=rf.get();
        std::cout<<"DONE "<<data.size()<<std::endl;
        
    }
    
    
    /*
     std::cout<<label<<std::endl;
     parsplice::Rd reft;
     std::shared_future<Rd> rf= ds.get( DDS_GET | DDS_FORWARD, 1, label);
     //ds.get( DDS_GET | DDS_FORWARD, 1, label);
     //ds.get( DDS_GET | DDS_FORWARD, 1, label);
     
     while( rf.wait_for(std::chrono::seconds(0)) != std::future_status::ready){
     
     }
     
     
     Rd data=rf.get();
     std::cout<<"DONE "<<data.size()<<std::endl;
     */
    /*
     for(int i=0;i<data.size();i++){
     std::cout<<int(data[i]);
     }
     */
    
    while(true){
        
    }
}
#endif



#if 0
////////////////////////////////////////////////////////////////////////////////
// test the whole setup
////////////////////////////////////////////////////////////////////////////////


//communicators
//intercomm with splicer and all dispatchers


int predictionTime=30;
double wMin=0.001;
MPI_Comm *depthComm;
MPI_Comm validatorComm;
MPI_Comm *dbComm;
MPI_Comm dispatchComm;
MPI_Comm dispatch2Comm;

//create a validator comm that is just a copy of world
MPI_Comm_dup(MPI_COMM_WORLD, &validatorComm);
//create a second dispatch comm that is just a copy of world
//the dispatchers and producers can talk over that comm
MPI_Comm_dup(MPI_COMM_WORLD, &dispatch2Comm);


if(depth==0){
    //splicer validator db
    
    //create the dispatchComm
    //dispatchComm includes the root and the layer. The broadcast of predictions it done over this comm.
    MPI_Comm_split(MPI_COMM_WORLD, 1, rank, &dispatchComm );
    
    ParSpliceSplicer splicer(dispatchComm, predictionTime, wMin);
    Validator validator(validatorComm,parent,children, 3, splicer.validatedDB,splicer.validatedLock);
    
    Visit v;
    v.label=1;
    splicer.officialTrajectory.appendVisit(v);
    splicer.loop();
}


if(depth==1){
    //validator db
    
    //create the dispatchComm
    MPI_Comm_split(MPI_COMM_WORLD, 1, rank, &dispatchComm );
    
    
    //MPI_Comm_split(MPI_COMM_WORLD, depth, rank, &depthComm);
    
    
    //these will remain unused here
    TrajDB validatedDB;
    pthread_mutex_t lock;
    pthread_mutex_init(&lock, NULL);
    Validator validator(validatorComm,parent,children, 3, validatedDB,lock);
    Dispatcher dispatcher(dispatchComm,dispatch2Comm, predictionTime, parent);
    
    
    while(true){
        
    }
}

if(depth==2){
    //producers
    
    boost::random::mt19937 rand;
    rand.seed(rank*1234);
    boost::random::uniform_int_distribution<int> d(1,10);
    std::vector<unsigned int> uintCommBuff(BUFFER_SIZE,0);
    
    
    //create the dispatchComm
    MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, rank, &dispatchComm );
    
    
    while(true){
        BOOST_LOG_SEV(lg,info) <<"#Request task from dispatcher";
        unsigned current;
        MPI_Send(&current, 1, MPI::UNSIGNED , parent, 1, dispatch2Comm);
        MPI_Status status;
        MPI_Recv(&current, 1, MPI::UNSIGNED , parent, MPI_ANY_TAG, dispatch2Comm, &status);
        
        BOOST_LOG_SEV(lg,info) <<"#Received task "<<current;
        
        TrajDB db;
        Trajectory t;
        Visit v;
        
        
        //generate a segment
        v.label=current;
        v.duration=d(rand);
        t.appendVisit(v);
        
        v.label=d(rand);
        v.duration=d(rand);
        v.duration=0;
        t.appendVisit(v);
        
        db[t.front().label].push_back(t);
        
        sleep(1);
        
        {
            logging::record rec = lg.open_record();
            if (rec)
            {
                logging::record_ostream strm(rec);
                strm<<"========== SENDING ==========\n";
                
                //loop over initial states
                for(TrajDB::iterator it=db.begin();it!=db.end();it++){
                    //loop over trajectories
                    for(std::deque<Trajectory>::iterator itt=it->second.begin();itt!=it->second.end();itt++){
                        strm<<"====================\n";
                        for(std::list<Visit>::iterator ittt=itt->visits.begin();ittt!=itt->visits.end();ittt++){
                            strm<<" -- ( "<<ittt->label<<" ) "<<ittt->duration;
                        }
                        strm<<"\n";
                    }
                }
                
                
                strm.flush();
                lg.push_record(boost::move(rec));
            }
            
            
        }
        
        int i=0;
        unsigned waste=d(rand);
        serializeDB(db, waste, uintCommBuff, i);
        MPI_Send(&(uintCommBuff[0]),i,MPI::UNSIGNED,parent,0,validatorComm);
        BOOST_LOG_SEV(lg,info) <<"#Sending reply "<<current;
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    while(true){
        
    }
    
}
#endif
////////////////////////////////////////////////////////////////////////////////
// test the Validator
////////////////////////////////////////////////////////////////////////////////
#if 0
if(rank==0){
    std::set<int> children;
    children.insert(1);
    //children.insert(2);
    TrajDB validated;
    
    Validator v(MPI_COMM_WORLD,-1,children, 3, validated);
    v.server();
}

if(rank==1){
    std::set<int> children;
    //children.insert(1);
    children.insert(2);
    children.insert(3);
    TrajDB validated;
    
    Validator v(MPI_COMM_WORLD,0,children, 3, validated);
    v.server();
}

if(rank>1){
    
    
    
    
    std::vector<unsigned int> uintCommBuff(BUFFER_SIZE,0);
    
    
    boost::random::mt19937 rand;
    rand.seed(rank*1234);
    
    
    boost::random::uniform_int_distribution<int> d(1,10);
    
    
    
    for(int l=0;l<20;l++){
        TrajDB out;
        for(int j=0;j<3;j++){
            Trajectory t;
            
            int n=d(rand);
            n=1;
            for(int k=0;k<n;k++){
                Visit v;
                int ind=d(rand);
                int duration=d(rand);
                v.label=ind;
                v.duration=duration;
                t.appendVisit(v);
            }
            out[t.front().label].push_back(t);
        }
        
        
        {
            logging::record rec = lg.open_record();
            if (rec)
            {
                logging::record_ostream strm(rec);
                strm<<"========== SENDING ==========\n";
                
                //loop over initial states
                for(TrajDB::iterator it=out.begin();it!=out.end();it++){
                    //loop over trajectories
                    for(std::deque<Trajectory>::iterator itt=it->second.begin();itt!=it->second.end();itt++){
                        strm<<"====================\n";
                        for(std::list<Visit>::iterator ittt=itt->visits.begin();ittt!=itt->visits.end();ittt++){
                            strm<<" -- ( "<<ittt->label<<" ) "<<ittt->duration;
                        }
                        strm<<"\n";
                    }
                }
                
                
                strm.flush();
                lg.push_record(boost::move(rec));
            }
            
            
        }
        
        //sleep(d(rand));
        
        int i=0;
        unsigned waste=10;
        serializeDB(out, waste, uintCommBuff, i);
        MPI_Send(&(uintCommBuff[0]),i,MPI::UNSIGNED,1,0,MPI_COMM_WORLD);
    }
    
    
    
    
    while(true){
        
    }
}
#endif

////////////////////////////////////////////////////////////////////////////////
// test the BKL
////////////////////////////////////////////////////////////////////////////////
#if 0
boost::random::mt19937 rand;
int minCount=1;
int maxCount=90;
int nStates=1000;
int order=10;

int nSample=1000000;

std::pair<Label, unsigned> out;
boost::random::uniform_int_distribution<int> d(minCount,maxCount);
boost::random::uniform_int_distribution<int> s(0,nStates-1);


//fill in the MarkovModel

for (int order=100;order<1000;order+=10){
    //std::cout<<"Bs"<<std::endl;
    MC mc;
    {
        //boost::timer::auto_cpu_timer t;
        
        for(unsigned int i=0;i<nStates;i++){
            
            for(int j=0;j<order;j++){
                unsigned int is=s(rand);
                int nt=d(rand);
                
                for(int k=0;k<nt;k++){
                    mc.incrementTrans(i,is);
                    mc.incrementTime(i,is*5);
                }
                
            }
        }
    }
    //std::cout<<"Be"<<std::endl;
    //std::cout<<"order: "<<order<<std::endl;
    
    
    {
        boost::timer::auto_cpu_timer t(3, "%w seconds\n");
        
        double w=1;
        
        Label cs=0;
        for(int i=0;i<nSample;i++){
            out=mc.sampleBKL(cs,w);
            
            cs=out.first;
            //std::cout<<cs<<" "<<out.second<<std::endl;
        }
    }
}
#endif

