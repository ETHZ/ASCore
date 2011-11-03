ROOTCFLAGS     = $(shell root-config --cflags)
ROOTLIBS       = $(shell root-config --libs)
ROOTGLIBS      = $(shell root-config --glibs)

INCLUDES       = -I./include -I$(CMSSW_RELEASE_BASE)/src/

CXX            = g++
CXXFLAGS       = -g -fPIC -Wno-deprecated -D_GNU_SOURCE -O2 $(INCLUDES) 
LD             = g++
LDFLAGS        = -g 
SOFLAGS        = -O --no_exceptions -shared


CXXFLAGS      += $(ROOTCFLAGS)
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) -lMinuit -lMinuit2 -lTreePlayer
GLIBS          = $(filter-out -lNew, $(NGLIBS)) -L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH) -lFWCoreFWLite -lCondFormatsJetMETObjects

SRCS           = src/base/TreeClassBase.C src/base/TreeReader.cc src/base/TreeAnalyzerBase.cc src/base/UserAnalysisBase.cc \
                 src/UserAnalyzer.cc src/TreeSkimmer.cc src/UserAnalysis.cc \
                 src/helper/PUWeight.C src/helper/AnaClass.cc src/helper/Davismt2.cc src/helper/LeptJetStat.cc src/helper/Hemisphere.cc src/helper/MetaTreeClassBase.C

OBJS           = $(patsubst %.C,%.o,$(SRCS:.cc=.o))

.SUFFIXES: .cc,.C,.hh,.h
.PHONY : clean purge all depend

# Rules ====================================
all: RunUserAnalyzer RunTreeSkimmer 

RunUserAnalyzer: src/exe/RunUserAnalyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

RunTreeSkimmer: src/exe/RunTreeSkimmer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

clean:
	$(RM) $(OBJS)	
	$(RM) RunUserAnalyzer
	$(RM) RunTreeSkimmer

purge:
	$(RM) $(OBJS)

deps: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it

src/base/TreeClassBase.o: ./include/base/TreeClassBase.h
src/base/TreeReader.o: ./include/base/TreeReader.hh
src/base/TreeReader.o: ./include/base/TreeClassBase.h
src/base/TreeAnalyzerBase.o: /usr/include/stdlib.h /usr/include/features.h
src/base/TreeAnalyzerBase.o: /usr/include/sys/cdefs.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/wordsize.h
src/base/TreeAnalyzerBase.o: /usr/include/gnu/stubs.h
src/base/TreeAnalyzerBase.o: /usr/include/gnu/stubs-64.h
src/base/TreeAnalyzerBase.o: /usr/include/sys/types.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/types.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/typesizes.h
src/base/TreeAnalyzerBase.o: /usr/include/time.h /usr/include/endian.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/endian.h
src/base/TreeAnalyzerBase.o: /usr/include/sys/select.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/select.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/sigset.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/time.h
src/base/TreeAnalyzerBase.o: /usr/include/sys/sysmacros.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/pthreadtypes.h
src/base/TreeAnalyzerBase.o: /usr/include/alloca.h
src/base/TreeAnalyzerBase.o: ./include/base/TreeAnalyzerBase.hh
src/base/TreeAnalyzerBase.o: ./include/base/TreeReader.hh
src/base/TreeAnalyzerBase.o: ./include/base/TreeClassBase.h
src/base/TreeAnalyzerBase.o: ./include/helper/Utilities.hh
src/base/TreeAnalyzerBase.o: /usr/include/stdio.h /usr/include/libio.h
src/base/TreeAnalyzerBase.o: /usr/include/_G_config.h /usr/include/wchar.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/stdio_lim.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/sys_errlist.h
src/base/TreeAnalyzerBase.o: /usr/include/assert.h
src/base/UserAnalysisBase.o: ./include/base/TreeReader.hh
src/base/UserAnalysisBase.o: ./include/base/TreeClassBase.h
src/base/UserAnalysisBase.o: /usr/include/stdlib.h /usr/include/features.h
src/base/UserAnalysisBase.o: /usr/include/sys/cdefs.h
src/base/UserAnalysisBase.o: /usr/include/bits/wordsize.h
src/base/UserAnalysisBase.o: /usr/include/gnu/stubs.h
src/base/UserAnalysisBase.o: /usr/include/gnu/stubs-64.h
src/base/UserAnalysisBase.o: /usr/include/sys/types.h
src/base/UserAnalysisBase.o: /usr/include/bits/types.h
src/base/UserAnalysisBase.o: /usr/include/bits/typesizes.h
src/base/UserAnalysisBase.o: /usr/include/time.h /usr/include/endian.h
src/base/UserAnalysisBase.o: /usr/include/bits/endian.h
src/base/UserAnalysisBase.o: /usr/include/sys/select.h
src/base/UserAnalysisBase.o: /usr/include/bits/select.h
src/base/UserAnalysisBase.o: /usr/include/bits/sigset.h
src/base/UserAnalysisBase.o: /usr/include/bits/time.h
src/base/UserAnalysisBase.o: /usr/include/sys/sysmacros.h
src/base/UserAnalysisBase.o: /usr/include/bits/pthreadtypes.h
src/base/UserAnalysisBase.o: /usr/include/alloca.h
src/base/UserAnalysisBase.o: ./include/helper/pdgparticle.hh
src/base/UserAnalysisBase.o: ./include/helper/Monitor.hh
src/base/UserAnalysisBase.o: ./include/base/UserAnalysisBase.hh
src/base/UserAnalysisBase.o: ./include/base/TreeReader.hh
src/base/UserAnalysisBase.o: ./include/helper/Utilities.hh
src/base/UserAnalysisBase.o: /usr/include/stdio.h /usr/include/libio.h
src/base/UserAnalysisBase.o: /usr/include/_G_config.h /usr/include/wchar.h
src/base/UserAnalysisBase.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/base/UserAnalysisBase.o: /usr/include/bits/stdio_lim.h
src/base/UserAnalysisBase.o: /usr/include/bits/sys_errlist.h
src/base/UserAnalysisBase.o: ./include/helper/PUWeight.h
src/base/UserAnalysisBase.o: /swshare/cms/slc5_amd64_gcc434/cms/cmssw/CMSSW_4_1_3/src/CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h
src/base/UserAnalysisBase.o: /swshare/cms/slc5_amd64_gcc434/cms/cmssw/CMSSW_4_1_3/src/CondFormats/JetMETObjects/interface/JetCorrectorParameters.h
src/base/UserAnalysisBase.o: /swshare/cms/slc5_amd64_gcc434/cms/cmssw/CMSSW_4_1_3/src/FWCore/Utilities/interface/Exception.h
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/is_base_and_derived.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/is_class.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/config.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/config.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/config/user.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/config/select_compiler_config.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/config/compiler/gcc.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/config/select_stdlib_config.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/config/select_platform_config.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/config/posix_features.hpp
src/base/UserAnalysisBase.o: /usr/include/unistd.h
src/base/UserAnalysisBase.o: /usr/include/bits/posix_opt.h
src/base/UserAnalysisBase.o: /usr/include/bits/confname.h
src/base/UserAnalysisBase.o: /usr/include/getopt.h
src/base/UserAnalysisBase.o: /usr/include/boost/config/suffix.hpp
src/base/UserAnalysisBase.o: /usr/include/limits.h
src/base/UserAnalysisBase.o: /usr/include/bits/posix1_lim.h
src/base/UserAnalysisBase.o: /usr/include/bits/local_lim.h
src/base/UserAnalysisBase.o: /usr/include/linux/limits.h
src/base/UserAnalysisBase.o: /usr/include/bits/posix2_lim.h
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/is_union.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/remove_cv.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/broken_compiler_spec.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/aux_/lambda_support.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/aux_/config/lambda.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/aux_/config/ttp.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/aux_/config/msvc.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/config.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/aux_/config/gcc.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/aux_/config/workaround.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/detail/workaround.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/aux_/config/ctps.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/detail/cv_traits_impl.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/detail/yes_no_type.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/detail/type_trait_def.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/detail/template_arity_spec.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/int.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/int_fwd.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/aux_/adl_barrier.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/aux_/config/adl.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/aux_/config/intel.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/aux_/nttp_decl.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/aux_/config/nttp.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/preprocessor/cat.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/preprocessor/config/config.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/aux_/integral_wrapper.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/integral_c_tag.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/aux_/config/static_constant.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/aux_/static_cast.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/aux_/template_arity_fwd.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/aux_/preprocessor/params.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/aux_/config/preprocessor.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/preprocessor/comma_if.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/preprocessor/punctuation/comma_if.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/preprocessor/control/if.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/preprocessor/control/iif.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/preprocessor/logical/bool.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/preprocessor/facilities/empty.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/preprocessor/punctuation/comma.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/preprocessor/repeat.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/preprocessor/repetition/repeat.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/preprocessor/debug/error.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/preprocessor/detail/auto_rec.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/preprocessor/tuple/eat.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/preprocessor/inc.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/preprocessor/arithmetic/inc.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/aux_/config/overload_resolution.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/detail/type_trait_undef.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/intrinsics.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/detail/bool_trait_def.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/integral_constant.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/bool.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/bool_fwd.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/integral_c.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/mpl/integral_c_fwd.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/detail/bool_trait_undef.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/detail/ice_and.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/detail/ice_not.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/is_same.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/is_convertible.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/is_array.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/add_reference.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/is_reference.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/ice.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/detail/ice_or.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/detail/ice_eq.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/is_arithmetic.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/is_integral.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/is_float.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/type_traits/is_abstract.hpp
src/base/UserAnalysisBase.o: /usr/include/boost/static_assert.hpp
src/UserAnalyzer.o: ./include/UserAnalyzer.hh
src/UserAnalyzer.o: ./include/base/TreeAnalyzerBase.hh
src/UserAnalyzer.o: ./include/base/TreeReader.hh
src/UserAnalyzer.o: ./include/base/TreeClassBase.h
src/UserAnalyzer.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/UserAnalyzer.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/UserAnalyzer.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/UserAnalyzer.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
src/UserAnalyzer.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/UserAnalyzer.o: /usr/include/_G_config.h /usr/include/wchar.h
src/UserAnalyzer.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/UserAnalyzer.o: /usr/include/bits/stdio_lim.h
src/UserAnalyzer.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/UserAnalyzer.o: /usr/include/sys/types.h /usr/include/time.h
src/UserAnalyzer.o: /usr/include/endian.h /usr/include/bits/endian.h
src/UserAnalyzer.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/UserAnalyzer.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/UserAnalyzer.o: /usr/include/sys/sysmacros.h
src/UserAnalyzer.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/UserAnalyzer.o: ./include/base/TreeReader.hh ./include/UserAnalysis.hh
src/UserAnalyzer.o: ./include/base/UserAnalysisBase.hh
src/UserAnalyzer.o: ./include/helper/pdgparticle.hh
src/UserAnalyzer.o: ./include/helper/PUWeight.h
src/UserAnalyzer.o: /swshare/cms/slc5_amd64_gcc434/cms/cmssw/CMSSW_4_1_3/src/CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h
src/TreeSkimmer.o: ./include/TreeSkimmer.hh
src/TreeSkimmer.o: ./include/base/TreeAnalyzerBase.hh
src/TreeSkimmer.o: ./include/base/TreeReader.hh
src/TreeSkimmer.o: ./include/base/TreeClassBase.h
src/TreeSkimmer.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/TreeSkimmer.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/TreeSkimmer.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/TreeSkimmer.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
src/TreeSkimmer.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/TreeSkimmer.o: /usr/include/_G_config.h /usr/include/wchar.h
src/TreeSkimmer.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/TreeSkimmer.o: /usr/include/bits/stdio_lim.h
src/TreeSkimmer.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/TreeSkimmer.o: /usr/include/sys/types.h /usr/include/time.h
src/TreeSkimmer.o: /usr/include/endian.h /usr/include/bits/endian.h
src/TreeSkimmer.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/TreeSkimmer.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/TreeSkimmer.o: /usr/include/sys/sysmacros.h
src/TreeSkimmer.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/TreeSkimmer.o: ./include/base/TreeReader.hh
src/UserAnalysis.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/UserAnalysis.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/UserAnalysis.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/UserAnalysis.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
src/UserAnalysis.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/UserAnalysis.o: /usr/include/_G_config.h /usr/include/wchar.h
src/UserAnalysis.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/UserAnalysis.o: /usr/include/bits/stdio_lim.h
src/UserAnalysis.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/UserAnalysis.o: /usr/include/sys/types.h /usr/include/time.h
src/UserAnalysis.o: /usr/include/endian.h /usr/include/bits/endian.h
src/UserAnalysis.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/UserAnalysis.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/UserAnalysis.o: /usr/include/sys/sysmacros.h
src/UserAnalysis.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/UserAnalysis.o: ./include/UserAnalysis.hh ./include/base/TreeReader.hh
src/UserAnalysis.o: ./include/base/TreeClassBase.h
src/UserAnalysis.o: ./include/base/UserAnalysisBase.hh
src/UserAnalysis.o: ./include/base/TreeReader.hh
src/UserAnalysis.o: ./include/helper/pdgparticle.hh
src/UserAnalysis.o: ./include/helper/PUWeight.h
src/UserAnalysis.o: /swshare/cms/slc5_amd64_gcc434/cms/cmssw/CMSSW_4_1_3/src/CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h
src/helper/PUWeight.o: ./include/helper/PUWeight.h
src/helper/AnaClass.o: ./include/helper/AnaClass.hh /usr/include/math.h
src/helper/AnaClass.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/helper/AnaClass.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/helper/AnaClass.o: /usr/include/gnu/stubs-64.h
src/helper/AnaClass.o: /usr/include/bits/huge_val.h
src/helper/AnaClass.o: /usr/include/bits/mathdef.h
src/helper/AnaClass.o: /usr/include/bits/mathcalls.h /usr/include/stdlib.h
src/helper/AnaClass.o: /usr/include/sys/types.h /usr/include/bits/types.h
src/helper/AnaClass.o: /usr/include/bits/typesizes.h /usr/include/time.h
src/helper/AnaClass.o: /usr/include/endian.h /usr/include/bits/endian.h
src/helper/AnaClass.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/helper/AnaClass.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/helper/AnaClass.o: /usr/include/sys/sysmacros.h
src/helper/AnaClass.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/helper/AnaClass.o: /usr/include/stdio.h /usr/include/libio.h
src/helper/AnaClass.o: /usr/include/_G_config.h /usr/include/wchar.h
src/helper/AnaClass.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/helper/AnaClass.o: /usr/include/bits/stdio_lim.h
src/helper/AnaClass.o: /usr/include/bits/sys_errlist.h
src/helper/AnaClass.o: ./include/helper/Utilities.hh
src/helper/AnaClass.o: ./include/helper/MetaTreeClassBase.h
src/helper/AnaClass.o: ./include/helper/Utilities.hh
src/helper/Davismt2.o: ./include/helper/Davismt2.h /usr/include/math.h
src/helper/Davismt2.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/helper/Davismt2.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/helper/Davismt2.o: /usr/include/gnu/stubs-64.h
src/helper/Davismt2.o: /usr/include/bits/huge_val.h
src/helper/Davismt2.o: /usr/include/bits/mathdef.h
src/helper/Davismt2.o: /usr/include/bits/mathcalls.h
src/helper/LeptJetStat.o: ./include/helper/LeptJetStat.h
src/helper/Hemisphere.o: ./include/helper/Hemisphere.hh
src/helper/Hemisphere.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/helper/Hemisphere.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/helper/Hemisphere.o: /usr/include/bits/wordsize.h
src/helper/Hemisphere.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
src/helper/Hemisphere.o: /usr/include/bits/types.h
src/helper/Hemisphere.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/helper/Hemisphere.o: /usr/include/_G_config.h /usr/include/wchar.h
src/helper/Hemisphere.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/helper/Hemisphere.o: /usr/include/bits/stdio_lim.h
src/helper/Hemisphere.o: /usr/include/bits/sys_errlist.h
src/helper/Hemisphere.o: /usr/include/stdlib.h /usr/include/sys/types.h
src/helper/Hemisphere.o: /usr/include/time.h /usr/include/endian.h
src/helper/Hemisphere.o: /usr/include/bits/endian.h /usr/include/sys/select.h
src/helper/Hemisphere.o: /usr/include/bits/select.h
src/helper/Hemisphere.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/helper/Hemisphere.o: /usr/include/sys/sysmacros.h
src/helper/Hemisphere.o: /usr/include/bits/pthreadtypes.h
src/helper/Hemisphere.o: /usr/include/alloca.h
src/helper/MetaTreeClassBase.o: ./include/helper/MetaTreeClassBase.h
