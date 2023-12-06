# echo "setup tautauAlg tautauAlg-00-00-01 in /besfs5/groups/psipp/psippgroup/public/xiaocong/bes3/7.0.3/workarea-7.0.3/Analysis"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmttautauAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmttautauAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=tautauAlg -version=tautauAlg-00-00-01 -path=/besfs5/groups/psipp/psippgroup/public/xiaocong/bes3/7.0.3/workarea-7.0.3/Analysis  -no_cleanup $* >${cmttautauAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=tautauAlg -version=tautauAlg-00-00-01 -path=/besfs5/groups/psipp/psippgroup/public/xiaocong/bes3/7.0.3/workarea-7.0.3/Analysis  -no_cleanup $* >${cmttautauAlgtempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmttautauAlgtempfile}
  unset cmttautauAlgtempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmttautauAlgtempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmttautauAlgtempfile}
unset cmttautauAlgtempfile
exit $cmtsetupstatus

