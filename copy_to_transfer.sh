#!/bin/bash
file=$1

# Some nice variables
api_key=`cat ~/transfer.api.key`
server="https://transfer.oicr.on.ca"

# Uploading the actual file and get attachment id for each file
attachment_id=`curl -X POST -H "Accept: application/json" --insecure --user "$api_key:x" -T "$file" https://transfer.oicr.on.ca/attachments/binary_upload?filename=${file} | jq .attachment.id`
echo $attachment_id

# Send the message
#cat <<EOF | curl -s -X POST --insecure --user "$api_key:x" -H "Content-Type: application/json" -d @- $server/link
#{"link":
#  {"attachment":"$attachment_id"}
#}
