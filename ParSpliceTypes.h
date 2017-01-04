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


#ifndef ParSplice_ParSpliceTypes_h
#define ParSplice_ParSpliceTypes_h

#include <stdio.h>
#include <vector>
#include <list>
#include <deque>
#include <memory>
#include <atomic>
#include <future>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/serialization/utility.hpp>
#include <boost/serialization/unordered_collections_save_imp.hpp>
#include <boost/serialization/unordered_collections_load_imp.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/serialization/deque.hpp>
#include <boost/serialization/list.hpp>



namespace parsplice{
    typedef unsigned Label;
    
    
    class AbstractSpinLock {
        
    public:
        virtual void lock()=0;
        virtual void unlock()=0;
    };
    
    class SpinLock : public AbstractSpinLock {
        std::atomic_flag locked = ATOMIC_FLAG_INIT ;
    public:
        void lock() {
            while (locked.test_and_set(std::memory_order_acquire)) { ; }
        }
        void unlock() {
            locked.clear(std::memory_order_release);
        }
    };
    
    class NoSpinLock : public AbstractSpinLock{
        
    public:
        void lock(){};
        void unlock(){};
    };
    
    
    
    
    struct Visit{
        friend class boost::serialization::access;
        Label label;
        unsigned int duration;
        
        Visit(){
            duration=0;
        };
        
        void serialize(std::vector<unsigned int> &v, int &i){
            v[i++]=(*this).label;
            v[i++]=(*this).duration;
        };
        
        void deserialize(std::vector<unsigned int> &v, int &i){
            (*this).label=v[i++];
            (*this).duration=v[i++];
        };
        
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){
            ar & label;
            ar & duration;
        };
        
    };
    
    
    struct Trajectory{
        friend class boost::serialization::access;
        
        
        Trajectory(){
            (*this).length=0;
        };
        
        
        void clear(){
            visits.clear();
            length=0;
            index=0;
        };
        
        //append a visit to a trajectory. Extends the last visit if possible
        void appendVisit(Visit &v, bool extend=true){
            (*this).length+=v.duration;
            
            if( extend and !(*this).empty() and (*this).back().label==v.label ){
                //extend the last visit
                (*this).back().duration+=v.duration;
                return;
            }
            //new visit
            (*this).visits.push_back(v);
        };
        
        
        //splice two trajectories
        bool splice(Trajectory &t){
            
            //do not splice empty trajectories
            if( (*this).empty() or t.empty() ){
                return false;
            }
            
            //splice only if the beginning and end match
            if( (*this).back().label == t.front().label ){
                (*this).back().duration+=t.front().duration;
                t.visits.pop_front();
                (*this).visits.splice((*this).visits.end(), t.visits);
                (*this).length+=t.length;
                
                //consume t
                t.visits.clear();
                t.length=0;
                return true;
            }
            
            //do nothing
            return false;
        };
        
        void serialize(std::vector<unsigned int> &v, int &i){
            v[i++]=(*this).index;
            v[i++]=(*this).visits.size();
            for(std::list<Visit>::iterator it=(*this).visits.begin();it!=(*this).visits.end();it++){
                it->serialize(v,i);
            }
        };
        
        void deserialize(std::vector<unsigned int> &v, int &i){
            (*this).length=0;
            (*this).visits.clear();
            (*this).index=v[i++];
            unsigned int nVisits=v[i++];
            for(int ii=0;ii<nVisits;ii++){
                Visit vis;
                vis.deserialize(v,i);
                (*this).length+=vis.duration;
                (*this).visits.push_back(vis);
            }
        };
        
        
        
        bool empty(){
            return visits.size()==0;
        };
        
        
        Visit& back(){
            return visits.back();
        };
        Visit& front(){
            return visits.front();
        };
        
        std::list<Visit> visits;
        unsigned int length;
        unsigned int index;
        
        
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){
            ar & visits;
            ar & length;
            ar & index;
        };
        
    };
    
    typedef boost::unordered_map< Label, std::deque<Trajectory> > TrajDB;
    
    typedef uint8_t Rdt;
    typedef std::vector<Rdt> Rd;
    
    
    class AbstractParSpliceWorkerDriver{
    public:
        virtual void generateSegment(parsplice::Label &label, Rd &ref, Rd &hot, parsplice::Trajectory &traj)=0;
        virtual void ref(parsplice::Label &label, parsplice::Rd &ref)=0;
    };
    
}


namespace boost {
    namespace serialization {
        
        template<class Archive, class Type, class Key, class Hash, class Compare, class Allocator >
        inline void save(
                         Archive & ar,
                         const boost::unordered_map<Key, Type, Hash, Compare, Allocator> &t,
                         const unsigned int /* file_version */
        ){
            boost::serialization::stl::save_unordered_collection(ar, t);
        }
        
        template<class Archive, class Type, class Key, class Hash, class Compare, class Allocator >
        inline void load(
                         Archive & ar,
                         boost::unordered_map<Key, Type, Hash, Compare, Allocator> &t,
                         const unsigned int /* file_version */
        ){
            
            boost::serialization::stl::load_unordered_collection<Archive,boost::unordered_map<Key, Type, Hash, Compare, Allocator>, boost::serialization::stl::archive_input_unordered_map<
            Archive,
            boost::unordered_map<Key, Type, Hash, Compare, Allocator>
            > >(ar,t);
            
    
        }
        
        // split non-intrusive serialization function member into separate
        // non intrusive save/load member functions
        template<class Archive, class Type, class Key, class Hash, class Compare, class Allocator >
        inline void serialize(
                              Archive & ar,
                              boost::unordered_map<Key, Type, Hash, Compare, Allocator> &t,
                              const unsigned int file_version
                              ){
            boost::serialization::split_free(ar, t, file_version);
        }
        
        
    } // serialization 
} // namespace boost 



#endif
