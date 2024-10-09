#docker pull --platform linux/amd64 nvidia/cuda:11.2.2-devel-ubuntu20.04
#docker run --rm -it --platform linux/amd64 nvidia/cuda:11.2.2-devel-ubuntu20.04 /bin/bash
FROM nvidia/cuda:11.2.2-devel-ubuntu20.04
LABEL authors="mainf"


ENTRYPOINT ["top", "-b"]