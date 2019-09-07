## `qserv`: managing Robo-AO's observing queue

Setup instructions:

```bash
docker build --rm -t qserv:latest -f Dockerfile .
docker run --name qserv -d --restart always -p 8081:8081 -v /path/to/Queue:/queue qserv:latest
#docker run -it --rm --name qserv -p 8081:8081 -v /Users/dmitryduev/web/qserv/_q:/queue qserv:latest
```

Make sure to replace `/path/to/Queue` with an actual system path. That's it!