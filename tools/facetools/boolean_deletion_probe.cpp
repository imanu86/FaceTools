// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

// Headless regression port of the original 1140-scenario scratch audit.
// Build/run only through run_boolean_deletion_probe.py; no GUI automation.
#include <catch2/catch.hpp>
#include "libslic3r/Model.hpp"
#include "libslic3r/FaceExtrusion.hpp"
#include "libslic3r/MeshBoolean.hpp"
#include "slic3r/GUI/GUI_ObjectList.hpp"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>
using namespace Slic3r;
using Slic3r::GUI::ItemForDelete;
using Slic3r::GUI::ItemType;
namespace {
bool fixed() { return std::getenv("MP_APPLY_SORT") && std::string(std::getenv("MP_APPLY_SORT"))=="1"; }
void apply_production_order(std::vector<ItemForDelete>& items) {
    if(fixed()) {
        // This include is extracted VERBATIM from the patched production TU.
        #include "production_delete_order.inc"
    }
}
bool execute_reverse(Model &model,const std::vector<ItemForDelete>& items,std::string &trace) {
    // Exact headless Model mutation performed by ObjectList's reverse loop.
    // Does not claim to execute wx callbacks/refresh/assembly/undo handling.
    for(auto it=items.rbegin();it!=items.rend();++it) {
        if(it->type!=ItemType::itVolume || it->obj_idx<0 || size_t(it->obj_idx)>=model.objects.size()) return false;
        ModelObject *o=model.objects[it->obj_idx];
        if(it->sub_obj_idx<0 || size_t(it->sub_obj_idx)>=o->volumes.size()) {trace+="OUT_OF_RANGE ";return false;}
        trace+=o->volumes[it->sub_obj_idx]->name+" ";
        o->delete_volume(size_t(it->sub_obj_idx));
    }
    return true;
}
std::vector<TriangleMesh> fixture() {
    TriangleMesh source;
    const char *fixture_path = std::getenv("FACETOOLS_PROBE_FIXTURE");
    REQUIRE(fixture_path != nullptr);
    REQUIRE(source.ReadSTLFile(fixture_path));
    std::vector<TriangleMesh> out{source};
    for(float z:{2.f,4.f}) {
        std::vector<int> ids;
        for(size_t i=0;i<source.its.indices.size();++i) {
            const Vec3i &f=source.its.indices[i];
            const Vec3f a=source.its.vertices[f[0]],b=source.its.vertices[f[1]],c=source.its.vertices[f[2]];
            if(a.z()==z&&b.z()==z&&c.z()==z&&(b-a).cross(c-a).z()>0) ids.push_back(int(i));
        }
        auto e=extrude_planar_face(source.its,ids,Vec3d(0,0,1.25));REQUIRE(e);out.emplace_back(std::move(e.mesh));
    }
    TriangleMesh fourth=make_cube(1,1,1);fourth.translate(Vec3f(-5,-5,0));out.push_back(std::move(fourth));
    return out;
}
TriangleMesh unite(const std::vector<TriangleMesh> &operands,size_t count) {
    TriangleMesh result=operands.front();
    for(size_t i=1;i<count;++i) {
        std::vector<TriangleMesh> out;int failed=0;
        MeshBoolean::mcut::make_boolean(result,operands[i],out,"UNION",nullptr,nullptr,[&]{++failed;});
        REQUIRE(failed==0);REQUIRE(out.size()==1);result=std::move(out.front());
    }
    REQUIRE(result.volume()==Approx(count==3?2625:2626));return result;
}
struct Snapshot {
    ModelVolume *ptr;ObjectID id;std::string name,material;int extruder;
    Transform3d transform;indexed_triangle_set mesh;
};
Snapshot capture(ModelVolume *v) {return {v,v->id(),v->name,v->material_id(),v->config.extruder(),v->get_matrix(),v->mesh().its};}
bool preserved(const ModelObject &o,const Snapshot &s,bool check_id=true) {
    if(std::find(o.volumes.begin(),o.volumes.end(),s.ptr)==o.volumes.end()) return false;
    ModelVolume *v=s.ptr;
    if(check_id && v->id()!=s.id) return false;
    if(v->name!=s.name || v->material_id()!=s.material || v->config.extruder()!=s.extruder || v->get_matrix().matrix()!=s.transform.matrix()) return false;
    const auto &m=v->mesh().its;
    if(m.vertices!=s.mesh.vertices || m.indices!=s.mesh.indices || m.properties.size()!=s.mesh.properties.size()) return false;
    return true;
}
}
TEST_CASE("Native deletion preserves real union across volume permutations and ID renewals", "[multipart_delete_audit]") {
    const auto meshes=fixture();size_t total=0,bad=0;
    for(size_t selected_count:{3u,4u}) {
        const TriangleMesh united=unite(meshes,selected_count);
        for(bool keep_unselected:{false,true}) {
            if(selected_count==4 && !keep_unselected) continue;
            std::vector<size_t> physical(selected_count+(keep_unselected?1:0));std::iota(physical.begin(),physical.end(),0);
            do {
                for(int renewal=0;renewal<3;++renewal) for(bool duplicates:{false,true}) {
                    Model model;REQUIRE_FALSE(model.is_need_backup());ModelObject *o=model.add_object();
                    std::vector<ModelVolume*> all,selected;
                    for(size_t i=0;i<selected_count;++i) {auto *v=o->add_volume(meshes[i],false);v->name=std::string(1,char('A'+i));all.push_back(v);selected.push_back(v);}
                    ModelVolume *unselected=nullptr;
                    if(keep_unselected) {
                        unselected=o->add_volume(make_cube(2,3,4),false);unselected->name="UNSELECTED";
                        unselected->config.set_key_value("extruder",new ConfigOptionInt(7));unselected->set_material_id("untouched_material");unselected->set_offset(Vec3d(23,17,8));
                        all.push_back(unselected);
                    }
                    for(size_t i=0;i<physical.size();++i)o->volumes[i]=all[physical[i]];
                    if(renewal==1)selected[1]->set_new_unique_id();
                    if(renewal==2){selected[0]->set_new_unique_id();selected.back()->set_new_unique_id();}
                    std::sort(selected.begin(),selected.end(),[](const ModelVolume *a,const ModelVolume *b){return a->id()<b->id();});
                    auto *result=o->add_volume(united,false);result->name="RESULT";result->config.set_key_value("extruder",new ConfigOptionInt(3));result->set_material_id("union_material");
                    const auto result_before=capture(result);const auto untouched_before=unselected?capture(unselected):result_before;
                    auto pos=std::find(o->volumes.begin(),o->volumes.end(),selected.front());REQUIRE(pos!=o->volumes.end());o->delete_volume(size_t(pos-o->volumes.begin()));
                    std::vector<ItemForDelete> items;
                    for(size_t i=1;i<selected.size();++i) {
                        pos=std::find(o->volumes.begin(),o->volumes.end(),selected[i]);REQUIRE(pos!=o->volumes.end());
                        items.emplace_back(ItemType::itVolume,0,int(pos-o->volumes.begin()));
                    }
                    if(duplicates)items.push_back(items.front());
                    apply_production_order(items);
                    std::string trace;const bool safe=execute_reverse(model,items,trace);
                    // ModelObject::delete_volume deliberately renews the final volume's ID.
                    // Pointer/geometry/metadata, not that renewed ID, identify the singleton.
                    const bool result_ok=preserved(*o,result_before,keep_unselected);
                    const bool unselected_ok=!keep_unselected||preserved(*o,untouched_before);
                    const bool ok=safe&&result_ok&&unselected_ok&&o->volumes.size()==(keep_unselected?2u:1u);
                    ++total;if(!ok)++bad;
                    std::ostringstream order;for(size_t i:physical)order<<i;
                    if(!ok && bad<=8)std::cout<<"DELETE_TRACE fixed="<<fixed()<<" order="<<order.str()<<" renewal="<<renewal<<" duplicates="<<duplicates<<" deletes="<<trace<<" result="<<result_ok<<" untouched="<<unselected_ok<<std::endl;
                    INFO("fixed="<<fixed()<<" order="<<order.str()<<" renewal="<<renewal<<" duplicates="<<duplicates<<" trace="<<trace);
                    CHECK(ok);
                }
            } while(std::next_permutation(physical.begin(),physical.end()));
        }
    }
    std::cout<<"DELETE_SUMMARY fixed="<<fixed()<<" cases="<<total<<" failed="<<bad<<std::endl;
}
TEST_CASE("Deletion items are volume-only and safe across objects and deferred orphans", "[multipart_delete_audit]") {
    std::vector<int> order{0,1,2,3,4};size_t total=0,bad=0;
    do {
        for(bool duplicates:{false,true}) {
            Model model;REQUIRE_FALSE(model.is_need_backup());std::vector<Snapshot> retained;
            std::vector<ItemForDelete> canonical;
            for(int oi=0;oi<3;++oi) {
                auto *o=model.add_object();
                const int count=oi==1?1:3;
                for(int vi=0;vi<count;++vi) {
                    auto *v=o->add_volume(make_cube(1+oi,1+vi,1),false);v->name="o"+std::to_string(oi)+"v"+std::to_string(vi);v->config.set_key_value("extruder",new ConfigOptionInt(2));
                    if(oi!=1&&vi==1)retained.push_back(capture(v));
                    else canonical.emplace_back(ItemType::itVolume,oi,vi);
                }
            }
            REQUIRE(canonical.size()==5);
            std::vector<ItemForDelete> items;for(int i:order)items.push_back(canonical[i]);if(duplicates)items.push_back(items.front());
            apply_production_order(items);std::string trace;bool ok=execute_reverse(model,items,trace);
            ok=ok&&model.objects.size()==3&&model.objects[1]->volumes.empty();
            ok=ok&&preserved(*model.objects[0],retained[0],false)&&preserved(*model.objects[2],retained[1],false);
            // Like ObjectList, orphan removal is deferred until all volume items finish.
            if(ok){model.delete_object(size_t(1));ok=model.objects.size()==2&&preserved(*model.objects[1],retained[1],false);}
            ++total;if(!ok)++bad;
            INFO("multi-object fixed="<<fixed()<<" duplicate="<<duplicates<<" trace="<<trace);CHECK(ok);
        }
    } while(std::next_permutation(order.begin(),order.end()));
    std::cout<<"MULTIOBJECT_SUMMARY fixed="<<fixed()<<" cases="<<total<<" failed="<<bad<<std::endl;
}
