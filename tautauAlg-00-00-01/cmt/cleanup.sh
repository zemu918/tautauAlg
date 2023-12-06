# echo "cleanup tautauAlg tautauAlg-00-00-01 in /besfs5/groups/psipp/psippgroup/public/xiaocong/bes3/7.0.3/workarea-7.0.3/Analysis"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmttautauAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmttautauAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=tautauAlg -version=tautauAlg-00-00-01 -path=/besfs5/groups/psipp/psippgroup/public/xiaocong/bes3/7.0.3/workarea-7.0.3/Analysis  $* >${cmttautauAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt cleanup -sh -pack=tautauAlg -version=tautauAlg-00-00-01 -path=/besfs5/groups/psipp/psippgroup/public/xiaocong/bes3/7.0.3/workarea-7.0.3/Analysis  $* >${cmttautauAlgtempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmttautauAlgtempfile}
  unset cmttautauAlgtempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmttautauAlgtempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmttautauAlgtempfile}
unset cmttautauAlgtempfile
return $cmtcleanupstatus

