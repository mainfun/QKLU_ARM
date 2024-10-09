#!/bin/bash

# ===========================
# 配置部分
# ===========================

# 本地要上传的文件路径
LOCAL_FILE="./bin/QKLU"

# 登录节点的用户和主机名/IP
LOGIN_USER="111.115.201.221"

# 登录节点上的目标目录
REMOTE_DIR="/home/wrq"

# 计算节点的用户和主机名/IP
COMPUTE_HOST="gpu01"

# 要在计算节点上执行的命令（这里假设上传的是一个脚本）
REMOTE_COMMAND="${REMOTE_DIR}/$(basename ${LOCAL_FILE})"

#执行参数
MTX_FILE="/home/wrq/qiankunLU/res/k3plates.mtx"

# ===========================
# 功能实现
# ===========================

# 1. 上传文件到登录节点
echo "上传文件到登录节点..."
scp "${LOCAL_FILE}" "${LOGIN_USER}:${REMOTE_DIR}/"

# 检查 scp 是否成功
# shellcheck disable=SC2181
if [ $? -ne 0 ]; then
    echo "文件上传失败。"
    exit 1
fi
echo "文件上传成功。"

# 2. 在登录节点上通过 SSH 登录到计算节点并执行文件
echo "在计算节点上执行文件..."
ssh ${LOGIN_USER} ssh ${COMPUTE_HOST} ${REMOTE_COMMAND} ${MTX_FILE}
echo "ssh ${LOGIN_USER} ssh ${COMPUTE_HOST} ${REMOTE_COMMAND} ${MTX_FILE}"

# 检查 SSH 命令是否成功
# shellcheck disable=SC2181
if [ $? -ne 0 ]; then
    echo "在计算节点上执行文件失败。"
    exit 1
fi
echo "文件在计算节点上成功执行。"

# 完成
echo "脚本执行完毕。"
